using SpinMC
using Test
using MPI
using Suppressor

#init MPI
MPI.Init()
comm = MPI.COMM_WORLD
commRank = MPI.Comm_rank(comm)
commSize = MPI.Comm_size(comm)
@test commSize == 2

#run tests
a1 = (3/2, sqrt(3)/2)
a2 = (3/2, -sqrt(3)/2)
uc = UnitCell(a1,a2)
b1 = addBasisSite!(uc, (0.0, 0.0))
addInteraction!(uc, b1, b1, [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0], (1, 0))
addInteraction!(uc, b1, b1, [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0], (0, 1))
addInteraction!(uc, b1, b1, [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0], (1, -1))
setField!(uc, b1, [0.0, 0.0, 0.5])
triangularlattice = Lattice(uc, (8,8))

if commRank == 0
    m_MPI = MonteCarlo(triangularlattice, 2.0, 10000, 1000000, seed=UInt(1))
    @suppress run!(m_MPI)
    e,e2 = means(m_MPI.observables.energy)
    @test isapprox(e, -2.9806, rtol=1e-3)
    @test isapprox(e2, 8.8882, rtol=1e-3)
else
    m_MPI = MonteCarlo(triangularlattice, 1.7, 10000, 1000000, seed=UInt(2))
    @suppress run!(m_MPI)
    e,e2 = means(m_MPI.observables.energy)
    @test isapprox(e, -2.8836, rtol=1e-3)
    @test isapprox(e2, 8.3209, rtol=1e-3)
end

#finalize MPI
MPI.Finalize()
@test MPI.Finalized()