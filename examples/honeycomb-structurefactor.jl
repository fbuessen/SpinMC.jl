using SpinMC
using LinearAlgebra

a1 = (3/2, sqrt(3)/2)
a2 = (3/2, -sqrt(3)/2)
uc = UnitCell(a1,a2)

b1 = addBasisSite!(uc, (0.0, 0.0))
b2 = addBasisSite!(uc, (1.0, 0.0))

M = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
addInteraction!(uc, b1, b2, M, (0, 0))
addInteraction!(uc, b1, b2, M, (0, -1))
addInteraction!(uc, b1, b2, M, (-1, 0))

L = (16, 16)
lattice = Lattice(uc, L)

thermalizationSweeps = 50000
measurementSweeps = 50000
beta = 10.0
m = MonteCarlo(lattice, beta, thermalizationSweeps, measurementSweeps)
run!(m)

# Fourier transform correlations to compute structure factor. 
N = 256
correlation = mean(m.observables.correlation) # The correlation is measured with respect to the first spin, i.e. the i-th entry is the correlation dot(S_1,S_i). 
kx = collect(range(-2pi,2pi,length=N))
ky = collect(range(-2pi,2pi,length=N))
structurefactor = zeros(N,N)
for i in 1:N
    for j in 1:N
        z = 0.0
        # Compute Fourier transformation at momentum (kx, ky). The real-space position of the i-th spin is obtained via getSitePosition(lattice,i). 
        for b in 1:length(lattice.unitcell.basis)
            for k in 1:length(lattice)
                z += cos(dot((kx[i],ky[j]),getSitePosition(lattice,k).-getSitePosition(lattice,b))) * correlation[k,b]
            end
        end
        structurefactor[j,i] = z / (length(lattice) * length(lattice.unitcell.basis))
    end
end

# Plot result
using Plots
heatmap(kx,ky,structurefactor,aspect_ratio=1,xrange=(-2pi,2pi),yrange=(-2pi,2pi),clims=(0,1),xlabel="kx", ylabel="ky")
xs = [4.0pi/3.0, 2.0pi/3.0, -2.0pi/3.0, -4.0pi/3.0, -2.0pi/3.0, 2.0pi/3.0, 4.0pi/3.0]
ys = [0.0, 2.0pi/sqrt(3.0), 2.0pi/sqrt(3.0), 0.0, -2.0pi/sqrt(3.0), -2.0pi/sqrt(3.0), 0.0]
plot!(xs, ys, label="Extended BZ")