using Suppressor

@testset "Integration tests" begin
    #cubic lattice antiferromagnet in a magnetic field
    a1 = (1.0, 0.0, 0.0)
    a2 = (0.0, 1.0, 0.0)
    a3 = (0.0, 0.0, 1.0)
    uc = UnitCell(a1,a2,a3)
    b1 = addBasisSite!(uc, (0.0, 0.0, 0.0))
    addInteraction!(uc, b1, b1, [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], (1, 0, 0))
    addInteraction!(uc, b1, b1, [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], (0, 1, 0))
    addInteraction!(uc, b1, b1, [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], (0, 0, 1))
    setField!(uc, b1, [0.0, 0.0, 0.2])
    cubiclattice = Lattice(uc, (4,4,4))

    m = MonteCarlo(cubiclattice, 2.0, 10000, 1000000, seed=UInt(0))
    @suppress run!(m)
    e,e2 = means(m.observables.energy)
    @test isapprox(e, -2.4850, rtol=1e-3)
    @test isapprox(e2, 6.1794, rtol=1e-3)

    #triangular lattice ferromagnet in a field
    a1 = (3/2, sqrt(3)/2)
    a2 = (3/2, -sqrt(3)/2)
    uc = UnitCell(a1,a2)
    b1 = addBasisSite!(uc, (0.0, 0.0))
    addInteraction!(uc, b1, b1, [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0], (1, 0))
    addInteraction!(uc, b1, b1, [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0], (0, 1))
    addInteraction!(uc, b1, b1, [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0], (1, -1))
    setField!(uc, b1, [0.0, 0.0, 0.5])
    triangularlattice = Lattice(uc, (8,8))

    m = MonteCarlo(triangularlattice, 2.0, 10000, 1000000, seed=UInt(0))
    @suppress run!(m)
    e,e2 = means(m.observables.energy)
    @test isapprox(e, -2.9806, rtol=1e-3)
    @test isapprox(e2, 8.8882, rtol=1e-3)

    #honeycomb lattice
    J1 = -1.1 #NN Heisenberg
    J2 = -0.3 #2NN Heisenberg
    K = -3.0 #NN Kitaev
    G = -0.5 #NN Gamma
    D = -0.6 #2NN Dzyaloshinskii-Moriya
    A = -0.8 #Onsite anisotropy

    M1x = [J1+K 0.0 0.0; 0.0 J1 G; 0.0 G J1]
    M1y = [J1 0.0 G; 0.0 J1+K 0.0; G 0.0 J1]
    M1z = [J1 G 0.0; G J1 0.0; 0.0 0.0 J1+K]
    M2 = [J2 D/sqrt(3.0) -D/sqrt(3.0); -D/sqrt(3.0) J2 D/sqrt(3.0); D/sqrt(3.0) -D/sqrt(3.0) J2]
    a1 = (3/2, sqrt(3)/2)
    a2 = (3/2, -sqrt(3)/2)
    uc = UnitCell(a1,a2)
    b1 = addBasisSite!(uc, (0.0, 0.0))
    b2 = addBasisSite!(uc, (1.0, 0.0))
    addInteraction!(uc, b1, b2, M1x, (0, 0))
    addInteraction!(uc, b1, b2, M1y, (0, -1))
    addInteraction!(uc, b1, b2, M1z, (-1, 0))
    addInteraction!(uc, b1, b1, M2, (1, 0))
    addInteraction!(uc, b1, b1, collect(transpose(M2)), (1, -1))
    addInteraction!(uc, b1, b1, M2, (0, -1))
    addInteraction!(uc, b2, b2, collect(transpose(M2)), (1, 0))
    addInteraction!(uc, b2, b2, M2, (1, -1))
    addInteraction!(uc, b2, b2, collect(transpose(M2)), (0, -1))
    setInteractionOnsite!(uc, b1, A*ones(3,3))
    setInteractionOnsite!(uc, b2, A*ones(3,3))
    honeycomblattice = Lattice(uc, (8,8))

    m = MonteCarlo(honeycomblattice, 0.5, 100000, 1000000, seed=UInt(0))
    @suppress run!(m)
    e,e2 = means(m.observables.energy)
    @test isapprox(e, -4.3855, rtol=1e-3)
    @test isapprox(e2, 19.2992, rtol=1e-3)
end