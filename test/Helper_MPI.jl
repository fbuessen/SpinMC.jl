using SpinMC
using Test
using MPI

#init MPI
MPI.Init()
comm = MPI.COMM_WORLD
commRank = MPI.Comm_rank(comm)
commSize = MPI.Comm_size(comm)
@test commSize == 2

#run tests
if commRank == 0
    rec = SpinMC.MPISendrecvFloat(42.0, 1, comm)
    @test rec == 43.0
elseif commRank == 1
    rec = SpinMC.MPISendrecvFloat(43.0, 0, comm)
    @test rec == 42.0
end

if commRank == 0
    SpinMC.MPISendBool(true, 1, comm)
    SpinMC.MPISendBool(false, 1, comm)
elseif commRank == 1
    rec = SpinMC.MPIRecvBool(0, comm)
    @test rec == true
    rec = SpinMC.MPIRecvBool(0, comm)
    @test rec == false
end

rec = SpinMC.MPIBcastBool(true, 0, comm)
@test rec == true
rec = SpinMC.MPIBcastBool(false, 0, comm)
@test rec == false

#finalize MPI
MPI.Finalize()
@test MPI.Finalized()