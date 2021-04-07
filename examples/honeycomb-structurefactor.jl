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

# Fourier transform correlations to compute structure factor
N = 256
correlation = mean(m.observables.correlation)
kx = collect(range(-2pi,2pi,length=N))
ky = collect(range(-2pi,2pi,length=N))
structurefactor = zeros(N,N)
for i in 1:N
    for j in 1:N
        z = 0.0
        for k in 1:length(lattice)
            z += cos(dot((kx[i],ky[j]),getSitePosition(lattice,k))) * correlation[k]
        end
        structurefactor[i,j] = z / length(lattice)
    end
end

# Plot result
using Plots
heatmap(kx,ky,structurefactor,aspect_ratio=1,xrange=(-2pi,2pi),yrange=(-2pi,2pi),clims=(0,1),xlabel="kx", ylabel="ky")
xs = [0.0, 2.0pi/sqrt(3.0), 2.0pi/sqrt(3.0), 0.0, -2.0pi/sqrt(3.0), -2.0pi/sqrt(3.0), 0.0]
ys = [4.0pi/3.0, 2.0pi/3.0, -2.0pi/3.0, -4.0pi/3.0, -2.0pi/3.0, 2.0pi/3.0, 4.0pi/3.0]
plot!(xs, ys, label="Extended BZ")