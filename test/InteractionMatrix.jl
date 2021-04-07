@testset "InteractionMatrix tests" begin
    @test_throws ErrorException SpinMC.InteractionMatrix(ones(4,3))

    m1 = SpinMC.InteractionMatrix()
    @test m1.m11 == 0.0
    @test m1.m12 == 0.0
    @test m1.m13 == 0.0
    @test m1.m21 == 0.0
    @test m1.m22 == 0.0
    @test m1.m23 == 0.0
    @test m1.m31 == 0.0
    @test m1.m32 == 0.0
    @test m1.m33 == 0.0

    M = rand(3,3)
    m2 = SpinMC.InteractionMatrix(M)
    @test [m2.m11 m2.m12 m2.m13; m2.m21 m2.m22 m2.m23; m2.m31 m2.m32 m2.m33] == M
end