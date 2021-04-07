@testset "Lattice tests" begin
    #test honeycomb lattice implementation
    a1 = (3/2, sqrt(3)/2)
    a2 = (3/2, -sqrt(3)/2)
    basis = [(0.0, 0.0), (1.0, 0.0)]
    uc = UnitCell(a1, a2)
    b1 = addBasisSite!(uc, basis[1])
    b2 = addBasisSite!(uc, basis[2])
    JKx = [2.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    JKy = [1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 1.0]
    JKz = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 2.0]
    addInteraction!(uc, b1, b2, JKx, (0, 0))
    addInteraction!(uc, b1, b2, JKy, (0, -1))
    addInteraction!(uc, b1, b2, JKz, (-1, 0))
    setInteractionOnsite!(uc, b1, 2.0*ones(3,3))
    setInteractionOnsite!(uc, b2, 3.0*ones(3,3))
    setField!(uc, b1, 4.0*ones(3))
    setField!(uc, b2, 5.0*ones(3))
    lattice = Lattice(uc, (4,4))

    @test typeof(lattice) == Lattice{2,3}
    @test size(lattice) == (4,4)
    @test length(lattice) == 4*4*2

    for i in 1:length(lattice) ; setSpin!(lattice, i, i.*(1.0,2.0,3.0)) ; end
    for i in 1:length(lattice) 
        @test getSpin(lattice, i) == i.*(1.0,2.0,3.0)
    end
    
    sitepositions = Vector{NTuple{2,Float64}}(undef,0)
    for n1 in 0:3
        for n2 in 0:3
            for b in 1:2
                push!(sitepositions, (n1 .* a1) .+ (n2 .* a2) .+ basis[b])
            end
        end
    end
    for i in 1:length(lattice)
        @test getSitePosition(lattice, i) in sitepositions

        #sites of basis b1
        getSitePosition(lattice, i) == (0.0, 0.0) && (@test SpinMC.getInteractionOnsite(lattice, i) == SpinMC.InteractionMatrix(2.0*ones(3,3)))
        getSitePosition(lattice, i) == (0.0, 0.0) && (@test SpinMC.getInteractionField(lattice, i) == Tuple(4.0*ones(3)))
        getSitePosition(lattice, i) == (4.0, 0.0) && (@test SpinMC.getInteractionOnsite(lattice, i) == SpinMC.InteractionMatrix(3.0*ones(3,3)))
        getSitePosition(lattice, i) == (4.0, 0.0) && (@test SpinMC.getInteractionField(lattice, i) == Tuple(5.0*ones(3)))

        #sites of basis b2
        getSitePosition(lattice, i) == (1.0, 0.0) && (@test SpinMC.getInteractionOnsite(lattice, i) == SpinMC.InteractionMatrix(3.0*ones(3,3)))
        getSitePosition(lattice, i) == (1.0, 0.0) && (@test SpinMC.getInteractionField(lattice, i) == Tuple(5.0*ones(3)))
        getSitePosition(lattice, i) == (3.0, 0.0) && (@test SpinMC.getInteractionOnsite(lattice, i) == SpinMC.InteractionMatrix(2.0*ones(3,3)))
        getSitePosition(lattice, i) == (3.0, 0.0) && (@test SpinMC.getInteractionField(lattice, i) == Tuple(4.0*ones(3)))

        interactionsites = SpinMC.getInteractionSites(lattice, i)
        interactionmatrices = SpinMC.getInteractionMatrices(lattice, i)
        
        for j in 1:length(interactionsites)
            offset = getSitePosition(lattice, interactionsites[j]) .- getSitePosition(lattice, i)
            if all(offset .≈ (1.0, 0.0))
                @test interactionmatrices[j] == SpinMC.InteractionMatrix(JKx)
            elseif all(offset .≈ (-1/2, sqrt(3)/2))
                @test interactionmatrices[j] == SpinMC.InteractionMatrix(JKy)
            elseif all(offset .≈ (-1/2, -sqrt(3)/2))
                @test interactionmatrices[j] == SpinMC.InteractionMatrix(JKz)
            end
        end
    end

    #test cubic lattice geometry
    a1 = (1.0, 0.0, 0.0)
    a2 = (0.0, 1.0, 0.0)
    a3 = (0.0, 0.0, 1.0)
    uc = UnitCell(a1, a2, a3)
    b = addBasisSite!(uc, (0.0,0.0,0.0))
    addInteraction!(uc, b, b, [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], (1,0,0))
    addInteraction!(uc, b, b, [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], (0,1,0))
    addInteraction!(uc, b, b, [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], (0,0,1))
    Lx = 2
    Ly = 3 
    Lz = 4
    lattice = Lattice(uc, (Lx, Ly, Lz))

    @test typeof(lattice) == Lattice{3,6}
    @test size(lattice) == (Lx, Ly, Lz)
    @test length(lattice) == Lx*Ly*Lz

    for i in 1:length(lattice)
        p = getSitePosition(lattice, i)
        neighbors = [p.+a1, p.-a1, p.+a2, p.-a2, p.+a3, p.-a3]
        for j in 1:6
            n = [ neighbors[j][k] for k in 1:3 ]
            n[1] < -0.5 && (n[1] += Lx)
            n[1] > Lx - 0.5 && (n[1] -= Lx)
            n[2] < -0.5 && (n[2] += Ly)
            n[2] > Ly - 0.5 && (n[2] -= Ly)
            n[3] < -0.5 && (n[3] += Lz)
            n[3] > Lz - 0.5 && (n[3] -= Lz)
            neighbors[j] = Tuple(n)
        end

        interactionsites = SpinMC.getInteractionSites(lattice, i)
        for j in 1:length(interactionsites)
            @test getSitePosition(lattice, interactionsites[j]) in neighbors
        end
    end
end