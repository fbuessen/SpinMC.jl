using Random
using LinearAlgebra

function uniformOnSphere(rng = Random.GLOBAL_RNG)
    phi = 2.0 * pi * rand(rng)
    z = 2.0 * rand(rng) - 1.0;
    r = sqrt(1.0 - z * z)
    return (r * cos(phi), r * sin(phi), z)
end

function exchangeEnergy(s1, M::InteractionMatrix, s2)::Float64
    return s1[1] * (M.m11 * s2[1] + M.m12 * s2[2] + M.m13 * s2[3]) + s1[2] * (M.m21 * s2[1] + M.m22 * s2[2] + M.m23 * s2[3]) + s1[3] * (M.m31 * s2[1] + M.m32 * s2[2] + M.m33 * s2[3])
end

function getEnergy(lattice::Lattice{D,N})::Float64 where {D,N}
    energy = 0.0

    for site in 1:length(lattice)
        s0 = getSpin(lattice, site)

        #two-spin interactions
        interactionSites = getInteractionSites(lattice, site)
        interactionMatrices = getInteractionMatrices(lattice, site)
        for i in 1:length(interactionSites)
            if site > interactionSites[i]
                energy += exchangeEnergy(s0, interactionMatrices[i], getSpin(lattice, interactionSites[i]))
            end
        end

        #onsite interaction
        energy += exchangeEnergy(s0, getInteractionOnsite(lattice, site), s0)

        #field interaction
        energy += dot(s0, getInteractionField(lattice, site))
    end

    return energy
end

function getEnergyDifference(lattice::Lattice{D,N}, site::Int, newState::Tuple{Float64,Float64,Float64})::Float64 where {D,N}
    dE = 0.0
    oldState = getSpin(lattice, site)
    ds = newState .- oldState

    #two-spin interactions
    interactionSites = getInteractionSites(lattice, site)
    interactionMatrices = getInteractionMatrices(lattice, site)
    for i in 1:length(interactionSites)
        dE += exchangeEnergy(ds, interactionMatrices[i], getSpin(lattice, interactionSites[i]))
    end

    #onsite interaction
    interactionOnsite = getInteractionOnsite(lattice, site)
    dE += exchangeEnergy(newState, interactionOnsite, newState) - exchangeEnergy(oldState, interactionOnsite, oldState)

    #field interaction
    dE += dot(ds, getInteractionField(lattice, site))

    return dE
end

function getMagnetization(lattice::Lattice{D,N}) where {D,N}
    mx, my, mz = 0.0, 0.0, 0.0
    for i in 1:length(lattice)
        spin = getSpin(lattice, i)
        mx += spin[1]
        my += spin[2]
        mz += spin[3]
    end
    return [mx, my, mz] / length(lattice)
end

function getCorrelation(lattice::Lattice{D,N}) where {D,N}
    corr = zeros(length(lattice), length(lattice.unitcell.basis))
    for i in 1:length(lattice.unitcell.basis)
        s0 = getSpin(lattice, i)
        for j in 1:length(lattice)
            corr[j,i] = dot(s0, getSpin(lattice, j))
        end
    end
    return corr
end