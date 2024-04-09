using Random
using Dates
using Printf
using MPI

mutable struct MonteCarloStatistics
    sweeps::Int
    attemptedLocalUpdates::Int
    acceptedLocalUpdates::Int
    attemptedReplicaExchanges::Int
    acceptedReplicaExchanges::Int
    initializationTime::Float64

    MonteCarloStatistics() = new(0, 0, 0, 0, 0, time())
end

mutable struct MonteCarlo{T<:Lattice,U<:AbstractRNG}
    lattice::T
    
    beta::Float64
    thermalizationSweeps::Int
    measurementSweeps::Int
    measurementRate::Int
    replicaExchangeRate::Int
    reportInterval::Int
    checkpointInterval::Int
    recursionIterations::Int
    recursionRate::Int

    rng::U
    seed::UInt
    sweep::Int
    recursion::Int

    observables::Observables
end

function MonteCarlo(
    lattice::T, 
    beta::Float64, 
    thermalizationSweeps::Int, 
    measurementSweeps::Int; 
    measurementRate::Int = 1, 
    replicaExchangeRate::Int = 10, 
    reportInterval::Int = round(Int, 0.05 * (thermalizationSweeps + measurementSweeps)), 
    checkpointInterval::Int = 3600, 
    recursionIterations::Int = 0, 
    recursionRate::Int = 1000, 
    rng::U = copy(Random.GLOBAL_RNG), 
    seed::UInt = rand(Random.RandomDevice(),UInt)
    ) where T<:Lattice where U<:AbstractRNG

    mc = MonteCarlo(deepcopy(lattice), beta, thermalizationSweeps, measurementSweeps, measurementRate, replicaExchangeRate, reportInterval, checkpointInterval, recursionIterations, recursionRate, rng, seed, 0, 0, Observables(lattice))
    Random.seed!(mc.rng, mc.seed)
    
    return mc
end

function run!(mc::MonteCarlo{T}; outfile::Union{String,Nothing}=nothing) where T<:Lattice
    #init MPI
    rank = 0
    commSize = 1
    allBetas = zeros(0)
    enableMPI = false
    if MPI.Initialized()
        commSize = MPI.Comm_size(MPI.COMM_WORLD)
        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        if commSize > 1
            allBetas = zeros(commSize)
            allBetas[rank + 1] = mc.beta
            MPI.Allgather!(UBuffer(allBetas, 1), MPI.COMM_WORLD)
            enableMPI = true
            rank == 0 && @printf("MPI detected. Enabling replica exchanges across %d simulations.\n", commSize)
        end

        #init parallell tempering recursion parameters
        if mc.recursionIterations > 0
            recursionStatistics = MonteCarloStatistics()
            weightedSum = 0
            weight = 0
        end
    end

    #init IO
    enableOutput = typeof(outfile) != Nothing
    if enableOutput
        enableMPI && (outfile *= "." * string(rank))
        isfile(outfile) && error("File ", outfile, " already exists. Terminating.")
    end
    
    #init spin configuration
    if mc.sweep == 0
        for i in 1:length(mc.lattice)
            setSpin!(mc.lattice, i, uniformOnSphere(mc.rng))
        end
    end

    #init Monte Carlo run
    totalSweeps = mc.thermalizationSweeps + mc.measurementSweeps
    partnerSpinConfiguration = deepcopy(mc.lattice.spins)
    energy = getEnergy(mc.lattice)

    #launch Monte Carlo run
    lastCheckpointTime = time()
    statistics = MonteCarloStatistics()
    rank == 0 && @printf("Simulation started on %s.\n\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))

    while mc.sweep < totalSweeps
        #perform local sweep
        for i in 1:length(mc.lattice)
            #select random spin
            site = rand(mc.rng, 1:length(mc.lattice))

            #propose new spin configuration
            newSpinState = uniformOnSphere(mc.rng)
            energyDifference = getEnergyDifference(mc.lattice, site, newSpinState)

            #check acceptance of new configuration
            statistics.attemptedLocalUpdates += 1
            p = exp(-mc.beta * energyDifference)
            if (rand(mc.rng) < min(1.0, p))
                setSpin!(mc.lattice, site, newSpinState)
                energy += energyDifference
                statistics.acceptedLocalUpdates += 1
            end
        end
        statistics.sweeps += 1

        #perform replica exchange
        if enableMPI && mc.sweep % mc.replicaExchangeRate == 0
            #determine MPI rank to exchagne configuration with
            if iseven(mc.sweep ÷ mc.replicaExchangeRate)
                partnerRank = iseven(rank) ? rank + 1 : rank - 1
            else
                partnerRank = iseven(rank) ? rank - 1 : rank + 1
            end

            if partnerRank >= 0 && partnerRank < commSize
                #obtain energy of new configuration
                partnerEnergy = MPISendrecvFloat(energy, partnerRank, MPI.COMM_WORLD)

                #check acceptance of new configuration
                statistics.attemptedReplicaExchanges += 1
                if mc.recursionIterations > 0
                    recursionStatistics.attemptedReplicaExchanges += 1
                end
                exchangeAccepted = false
                if iseven(rank)
                    p = exp(-(allBetas[rank + 1] - allBetas[partnerRank + 1]) * (partnerEnergy - energy))
                    exchangeAccepted = (rand(mc.rng) < min(1.0, p)) ? true : false
                    MPISendBool(exchangeAccepted, partnerRank, MPI.COMM_WORLD)
                else
                    exchangeAccepted = MPIRecvBool(partnerRank, MPI.COMM_WORLD)
                end
                if (exchangeAccepted)
                    energy = partnerEnergy
                    MPI.Sendrecv!(mc.lattice.spins, partnerRank, 0, partnerSpinConfiguration, partnerRank, 0, MPI.COMM_WORLD)
                    (mc.lattice.spins, partnerSpinConfiguration) = (partnerSpinConfiguration, mc.lattice.spins)
                    statistics.acceptedReplicaExchanges += 1
                    if mc.recursionIterations > 0
                        recursionStatistics.acceptedReplicaExchanges += 1
                    end
                end
            end
        end

        #perform measurement
        if mc.sweep >= mc.thermalizationSweeps
            if mc.sweep % mc.measurementRate == 0
                performMeasurements!(mc.observables, mc.lattice, energy)
            end
        end

        #parallel tempering temperature recursion
        if enableMPI && mc.sweep < mc.thermalizationSweeps && mc.recursionIterations > 0 && mc.recursion <= mc.recursionIterations && (mc.sweep + 1) % mc.recursionRate == 0
            #measure exchange rates
            replicaExchangeAcceptanceRate = recursionStatistics.acceptedReplicaExchanges / recursionStatistics.attemptedReplicaExchanges
            allReplicaExchangeAcceptanceRate = zeros(commSize)
            allReplicaExchangeAcceptanceRate[rank + 1] = replicaExchangeAcceptanceRate

            #share the exchange rates among processes
            MPI.Allgather!(UBuffer(allReplicaExchangeAcceptanceRate, 1), MPI.COMM_WORLD)

            #save info for recursion
            allReplicaExchangeAcceptanceRateMin = minimum(allReplicaExchangeAcceptanceRate)
            weightedSum += allReplicaExchangeAcceptanceRateMin * allBetas[rank + 1]
            weight += allReplicaExchangeAcceptanceRateMin

            if mc.recursion != mc.recursionIterations
                if rank == 0
                    #print statistics
                    str = ""
                    str *= @sprintf("Temperature readjustment %d / %d\n", mc.recursion + 1, mc.recursionIterations)
                    for n in 1:commSize
                        str *= @sprintf("\t\tsimulation %d replica exchange acceptance rate : %.2f%%\n", n - 1, 100.0 * allReplicaExchangeAcceptanceRate[n])
                    end
                    str *= @sprintf("\n")
                    print(str)

                    #calculate the new betas in root process
                    for n in 1:commSize
                        if allReplicaExchangeAcceptanceRate[n] == 0.0
                            allReplicaExchangeAcceptanceRate[n] = 0.01 # minimum allowed acceptance
                        end
                    end

                    #calculate auxiliary quantities
                    lambdaDenominator = 0.0
                    for n in 2:commSize
                        lambdaDenominator += allReplicaExchangeAcceptanceRate[n]*(allBetas[n] - allBetas[n - 1])
                    end
                    lambda = (allBetas[end]-allBetas[1]) / lambdaDenominator

                    #apply recursion
                    allBetasNew = zeros(commSize)
                    allBetasNew[1] = allBetas[1]
                    for n in 2:commSize
                        allBetasNew[n] = allBetasNew[n - 1] + lambda * allReplicaExchangeAcceptanceRate[n] * (allBetas[n] - allBetas[n - 1])
                    end
                    allBetas = copy(allBetasNew)
                end
                #share betas with rest of processes
                MPI.Bcast!(allBetas, 0, MPI.COMM_WORLD)

                #update betas
                mc.beta = allBetas[rank + 1]
            else
                #check if in every iteration the minimum acceptance was zero
                if weight == 0.0
                    error("Recursion failed, at least one replica exchange acceptance ratio was zero in every recursion iteration.")
                else
                    rank == 0 && @printf("\t\tRecursion successful\n\n")
                    #combine previous estimates of beta
                    allBetas[rank + 1] = weightedSum / weight
                    #update beta accross processes
                    MPI.Allgather!(UBuffer(allBetas, 1), MPI.COMM_WORLD)
                    mc.beta = allBetas[rank + 1]
                end
            end

            #reset statistics
            recursionStatistics = MonteCarloStatistics()

            mc.recursion += 1
        end

        #increment sweep
        statistics.sweeps += 1
        mc.sweep += 1

        #runtime statistics
        t = time()
        if mc.sweep % mc.reportInterval == 0
            #collect statistics
            progress = 100.0 * mc.sweep / totalSweeps
            thermalized = (mc.sweep >= mc.thermalizationSweeps) ? "YES" : "NO"
            sweeprate = statistics.sweeps / (t - statistics.initializationTime)
            sweeptime = 1.0 / sweeprate
            eta = (totalSweeps - mc.sweep) / sweeprate

            localUpdateAcceptanceRate = 100.0 * statistics.acceptedLocalUpdates / statistics.attemptedLocalUpdates
            if enableMPI
                replicaExchangeAcceptanceRate = 100.0 * statistics.acceptedReplicaExchanges / statistics.attemptedReplicaExchanges
                allLocalAppectanceRate = zeros(commSize)
                allLocalAppectanceRate[rank + 1] = localUpdateAcceptanceRate
                MPI.Allgather!(UBuffer(allLocalAppectanceRate, 1), MPI.COMM_WORLD)
                allReplicaExchangeAcceptanceRate = zeros(commSize)
                allReplicaExchangeAcceptanceRate[rank + 1] = replicaExchangeAcceptanceRate
                MPI.Allgather!(UBuffer(allReplicaExchangeAcceptanceRate, 1), MPI.COMM_WORLD)
            end

            #print statistics
            if rank == 0
                str = ""
                str *= @sprintf("Sweep %d / %d (%.1f%%)", mc.sweep, totalSweeps, progress)
                str *= @sprintf("\t\tETA : %s\n", Dates.format(Dates.now() + Dates.Second(round(Int64,eta)), "dd u yyyy HH:MM:SS"))
                str *= @sprintf("\t\tthermalized : %s\n", thermalized)
                str *= @sprintf("\t\tsweep rate : %.1f sweeps/s\n", sweeprate)
                str *= @sprintf("\t\tsweep duration : %.3f ms\n", sweeptime * 1000)
                
                if enableMPI
                    for n in 1:commSize
                        str *= @sprintf("\t\tsimulation %d update acceptance rate: %.2f%%\n", n - 1, allLocalAppectanceRate[n])
                        str *= @sprintf("\t\tsimulation %d replica exchange acceptance rate : %.2f%%\n", n - 1, allReplicaExchangeAcceptanceRate[n])
                    end
                else
                    str *= @sprintf("\t\tupdate acceptance rate: %.2f%%\n", localUpdateAcceptanceRate)
                end
                str *= @sprintf("\n")
                print(str)
            end

            #reset statistics
            statistics = MonteCarloStatistics()
        end

        #write checkpoint
        if enableOutput
            checkpointPending = time() - lastCheckpointTime >= mc.checkpointInterval
            enableMPI && (checkpointPending = MPIBcastBool(checkpointPending, 0, MPI.COMM_WORLD))
            if checkpointPending
                writeMonteCarlo(outfile, mc)
                lastCheckpointTime = time()
                rank == 0 && @printf("Checkpoint written on %s.\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
            end
        end
    end

    #write final checkpoint
    if enableOutput
        writeMonteCarlo(outfile, mc)
        rank == 0 && @printf("Checkpoint written on %s.\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
    end
    
    #return
    rank == 0 && @printf("Simulation finished on %s.\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
    return nothing    
end