using Suppressor

@testset "MonteCarlo tests" begin
    a1 = (3/2, sqrt(3)/2)
    a2 = (3/2, -sqrt(3)/2)
    basis = [(0.0, 0.0), (1.0, 0.0)]
    uc = UnitCell(a1, a2)
    b1 = addBasisSite!(uc, basis[1])
    b2 = addBasisSite!(uc, basis[2])
    addInteraction!(uc, b1, b2, [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0], (0, 0))
    addInteraction!(uc, b1, b2, [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0], (0, -1))
    addInteraction!(uc, b1, b2, [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0], (-1, 0))
    lattice = Lattice(uc, (4,4))

    beta = 1.0
    thermalizationSweeps = 100
    measurementSweeps = 200
    measurementRate = 2
    replicaExchangeRate = 3
    reportInterval = 50
    checkpointInterval = 1800
    seed = UInt(42)
    m1 = MonteCarlo(lattice, beta, thermalizationSweeps, measurementSweeps, measurementRate=measurementRate, replicaExchangeRate=replicaExchangeRate, reportInterval=reportInterval, checkpointInterval=checkpointInterval, seed=seed)
    m2 = MonteCarlo(lattice, beta, thermalizationSweeps, measurementSweeps, measurementRate=measurementRate, replicaExchangeRate=replicaExchangeRate, reportInterval=reportInterval, checkpointInterval=checkpointInterval, seed=seed)

    @test m1.beta == beta
    @test m1.thermalizationSweeps == thermalizationSweeps
    @test m1.measurementSweeps == measurementSweeps
    @test m1.measurementRate == measurementRate
    @test m1.replicaExchangeRate == replicaExchangeRate
    @test m1.reportInterval == reportInterval
    @test m1.checkpointInterval == checkpointInterval
    @test m1.seed == seed
    @test m1.sweep == 0

    @suppress run!(m1)
    @suppress run!(m2)
    @test m1.lattice.spins == m2.lattice.spins
    @test length(m1.observables.energy) == measurementSweeps / measurementRate
    @test length(m1.observables.magnetization) == measurementSweeps / measurementRate
    @test length(m1.observables.correlation) == measurementSweeps / measurementRate

    m3 = MonteCarlo(lattice, beta, thermalizationSweeps, measurementSweeps - 50, measurementRate=measurementRate, replicaExchangeRate=replicaExchangeRate, reportInterval=reportInterval, checkpointInterval=checkpointInterval, seed=seed)
    @suppress run!(m3)
    m3.measurementSweeps = measurementSweeps
    @suppress run!(m3)
    @test m1.lattice.spins == m3.lattice.spins
end