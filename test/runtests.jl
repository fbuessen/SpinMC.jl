using SpinMC
using Test

include("UnitCell.jl")
include("InteractionMatrix.jl")
include("Lattice.jl")
include("Observables.jl")
include("Spin.jl")
include("MonteCarlo.jl")
include("IntegrationTest.jl")

using MPI
nprocs = 2
testdir = @__DIR__
@testset "$f" for f in ["Helper_MPI.jl", "IntegrationTest_MPI.jl"]
    mpiexec() do cmd
        run(`$cmd -n $nprocs $(Base.julia_cmd()) $(joinpath(testdir, f))`)
        @test true
    end
end