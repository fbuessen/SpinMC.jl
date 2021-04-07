using SpinMC

#set up MPI
using MPI
MPI.Initialized() || MPI.Init()
commSize = MPI.Comm_size(MPI.COMM_WORLD)
commRank = MPI.Comm_rank(MPI.COMM_WORLD)
commSize > 1 || error("This example should be run on an MPI communicator of two ranks or more.")

#build lattice
a1 = (3/2, sqrt(3)/2)
a2 = (3/2, -sqrt(3)/2)
uc = UnitCell(a1,a2)

b1 = addBasisSite!(uc, (0.0, 0.0))
b2 = addBasisSite!(uc, (1.0, 0.0))

M = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
addInteraction!(uc, b1, b2, M, (0, 0))
addInteraction!(uc, b1, b2, M, (0, -1))
addInteraction!(uc, b1, b2, M, (-1, 0))

L = (16,16)
lattice = Lattice(uc, L)

#init MC simulation
thermalizationSweeps = 50000
measurementSweeps = 50000

tmin = 5.0
tmax = 10.0
beta = (commSize == 1) ? 1.0/tmin : 1.0 / (reverse([ tmax * (tmin / tmax)^(n/(commSize-1)) for n in 0:commSize-1 ])[commRank+1])

m = MonteCarlo(lattice, beta, thermalizationSweeps, measurementSweeps)

#run simulation
run!(m, outfile="simulation.h5")