@testset "UnitCell tests" begin
    @test typeof(UnitCell((1.0,))) == UnitCell{1}
    @test typeof(UnitCell((1.0,0.0),(0.0,1.0))) == UnitCell{2}
    @test typeof(UnitCell((1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0))) == UnitCell{3}
    @test typeof(UnitCell((1.0,0.0,0.0,0.0),(0.0,1.0,0.0,0.0),(0.0,0.0,1.0,0.0),(0.0,0.0,0.0,1.0))) == UnitCell{4}

    uc = UnitCell((1.0,0.0),(0.0,1.0))
    b1 = addBasisSite!(uc, (0.0, 0.0))
    b2 = addBasisSite!(uc, (0.5, 0.5))
    @test length(uc.basis) == 2
    @test length(uc.interactionsOnsite) == 2
    @test length(uc.interactionsField) == 2

    @test_throws ErrorException addInteraction!(uc, b1, b1, ones(3,3), (0,0))
    @test_throws ErrorException addInteraction!(uc, b1, b2, ones(4,3))
    addInteraction!(uc, b1, b1, ones(3,3), (1,0))
    addInteraction!(uc, b1, b1, ones(3,3), (0,1))
    @test length(uc.interactions) == 2

    @test_throws ErrorException setInteractionOnsite!(uc, b1, ones(4,3))
    setInteractionOnsite!(uc, b1, ones(3,3))
    setInteractionOnsite!(uc, b2, 2.0*ones(3,3))
    @test uc.interactionsOnsite[b1] == ones(3,3)
    @test uc.interactionsOnsite[b2] == 2.0*ones(3,3)

    @test_throws ErrorException setField!(uc, b1, ones(4))
    setField!(uc, b1, ones(3))
    setField!(uc, b2, 2.0*ones(3))
    @test uc.interactionsField[b1] == ones(3)
    @test uc.interactionsField[b2] == 2.0*ones(3)
end