using BinningAnalysis

mutable struct Observables
    energy::ErrorPropagator{Float64,32}
    magnetization::LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    magnetizationVector::LogBinner{Vector{Float64},32,BinningAnalysis.Variance{Vector{Float64}}}
    correlation::LogBinner{Array{Float64,2},32,BinningAnalysis.Variance{Array{Float64,2}}}
end

function Observables(lattice::T) where T<:Lattice
    return Observables(ErrorPropagator(Float64), LogBinner(Float64), LogBinner(zeros(Float64,3)), LogBinner(zeros(Float64,lattice.length, length(lattice.unitcell.basis)))) 
end

function performMeasurements!(observables::Observables, lattice::T, energy::Float64) where T<:Lattice
    #measure energy and energy^2
    push!(observables.energy, energy / length(lattice), energy * energy / (length(lattice) * length(lattice)))

    #measure magnetization
    m = getMagnetization(lattice)
    push!(observables.magnetization, norm(m))
    push!(observables.magnetizationVector, m)

    #measure spin correlations
    push!(observables.correlation, getCorrelation(lattice))
end