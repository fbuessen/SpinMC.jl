@testset "Observables tests" begin
    a1 = (3/2, sqrt(3)/2)
    a2 = (3/2, -sqrt(3)/2)
    basis = [(0.0, 0.0), (1.0, 0.0)]
    uc = UnitCell(a1, a2)
    b1 = addBasisSite!(uc, basis[1])
    b2 = addBasisSite!(uc, basis[2])
    addInteraction!(uc, b1, b2, [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], (0, 0))
    addInteraction!(uc, b1, b2, [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], (0, -1))
    addInteraction!(uc, b1, b2, [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], (-1, 0))
    setField!(uc, b1, [0.0, 0.0, -1.0])
    setField!(uc, b2, [0.0, 0.0, -1.0])
    lattice = Lattice(uc, (4,4))
    for i in 1:length(lattice)
        setSpin!(lattice, i, (0.0,0.0,1.0))
    end

    obs = Observables(lattice)
    energy = -42.0

    SpinMC.performMeasurements!(obs, lattice, energy)

    @test means(obs.energy)[1] ≈ energy / length(lattice)
    @test means(obs.energy)[2] ≈ (energy / length(lattice))^2
    @test mean(obs.magnetization) ≈ 1.0
    @test mean(obs.magnetizationVector) ≈ [0.0, 0.0, 1.0]
    @test all(mean(obs.correlation) .≈ ones(length(lattice)))
end