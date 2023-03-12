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

function Base.:*(M::SpinMC.InteractionMatrix, s)
    return ((M.m11 * s[1] + M.m12 * s[2] + M.m13 * s[3]) , (M.m21 * s[1] + M.m22 * s[2] + M.m23 * s[3]) , (M.m31 * s[1] + M.m32 * s[2] + M.m33 * s[3]))
end