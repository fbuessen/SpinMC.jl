@testset "Spin tests" begin
    x = (0.0, 0.0, 0.0)
    N = 1000000
    for i in 1:N
        x = x .+ SpinMC.uniformOnSphere()
    end
    @test all(isapprox.(x ./ N, (0.0, 0.0, 0.0), atol=1.0e-2))

    M = [1.1 2.2 3.3; 4.4 5.5 6.6; 7.7 8.8 9.9]
    s1 = (1.0, 2.0, 3.0)
    s2 = (4.0, 5.0, 6.0)
    @test SpinMC.exchangeEnergy(s1, SpinMC.InteractionMatrix(M), s2) ≈ 607.2

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
    for i in 1:length(lattice)
        setSpin!(lattice, i, (0.0,0.0,1.0))
    end
    @test getEnergy(lattice) ≈ -48.0
    @test SpinMC.getEnergyDifference(lattice, 1, (0.0, 0.0, -1.0)) ≈ 6.0
    setSpin!(lattice, 1, (0.0, 0.0, -1.0))
    @test getEnergy(lattice) ≈ -42.0

    @test getMagnetization(lattice) ≈ [0.0, 0.0, 30.0/32.0]
    @test all(getCorrelation(lattice)[:,1] .≈ [(i == 1 ? 1.0 : -1.0) for i in 1:32])
    @test all(getCorrelation(lattice)[:,2] .≈ [(i == 1 ? -1.0 : 1.0) for i in 1:32])
end