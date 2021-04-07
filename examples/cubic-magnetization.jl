using SpinMC

a1 = (1.0, 0.0, 0.0)
a2 = (0.0, 1.0, 0.0)
a3 = (0.0, 0.0, 1.0)
uc = UnitCell(a1,a2,a3)

b = addBasisSite!(uc, (0.0, 0.0, 0.0))
M = [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]
addInteraction!(uc, b, b, M, (1, 0, 0))
addInteraction!(uc, b, b, M, (0, 1, 0))
addInteraction!(uc, b, b, M, (0, 0, 1))

L = (8, 8, 8)
lattice = Lattice(uc, L)
thermalizationSweeps = 10000
measurementSweeps = 100000

tmin = 0.1
tmax = 10.0
N = 32
temperature = zeros(N)
heat = zeros(N)
dheat = zeros(N)
magnetization = zeros(N)
dmagnetization = zeros(N)
for i in 1:N
    println("Performing MC simulation ", i, "/", N)
    beta = 1.0 / (tmax * (tmin / tmax)^((i-1)/(N-1)))
    m = MonteCarlo(lattice, beta, thermalizationSweeps, measurementSweeps, reportInterval=20000)
    run!(m)
    temperature[i] = 1.0 / beta

    # Calculate specific heat 
    # (Note that if file output is activated in the run!() call, the specific heat is automatically included.)
    c(e) = beta * beta * (e[2] - e[1] * e[1]) * length(m.lattice)
    ∇c(e) = [-2.0 * beta * beta * e[1] * length(m.lattice), beta * beta * length(m.lattice)]
    heat[i] = mean(m.observables.energy, c)
    dheat[i] = std_error(m.observables.energy, ∇c)
    # Calculate magnetization
    magnetization[i] = mean(m.observables.magnetization)
    dmagnetization[i] = std_error(m.observables.magnetization)
end

using Plots
display(plot(temperature, heat, yerror=dheat, xaxis=:log, yrange=(0.0,2.5), xrange=(0.1,10.0), xlabel="temperature", ylabel="specific heat", label="8*8*8 cubic lattice"))
display(plot(temperature, magnetization, yerror=dmagnetization, xaxis=:log, yrange=(0.0,1.0), xrange=(0.1,10.0), xlabel="temperature", ylabel="magnetization", label="8*8*8 cubic lattice"))