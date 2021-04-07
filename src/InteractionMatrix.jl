struct InteractionMatrix
    m11::Float64
    m12::Float64
    m13::Float64
    m21::Float64
    m22::Float64
    m23::Float64
    m31::Float64
    m32::Float64
    m33::Float64
end

function InteractionMatrix()
    return InteractionMatrix(zeros(Float64,9)...)
end

function InteractionMatrix(M::T) where T<:AbstractMatrix
    size(M) == (3,3) || error("Interaction matrix must be of size 3x3")
    m = InteractionMatrix(M[1,1], M[1,2], M[1,3], M[2,1], M[2,2], M[2,3], M[3,1], M[3,2], M[3,3])
    return m
end