
"""
    acc_proj_power_method(A::SparseMatrixCSC{Int64,Int64}, X::Array{Float64,2}, n_it::Int64, t::Float64)

Runs the accelarated projected power iteration on the modularity matrix.

# Arguments
- `A::SparseMatrixCSC{Int64,Int64}` sparse graph adjacency matrix.
- `X::Array{Float64,2}` initial matrix.
- `n_it::Int64` maximum number of iterations.
- `t::Float64` tolerance on relative variation of consecutive objective values.
# Output
- `X::Array{Float64,2}: approximate stationary matrix of the iteration.

"""
function acc_proj_power_method(A::SparseMatrixCSC{Int64,Int64}, X::Array{Float64,2}, n_it::Int64, t::Float64)

    d::Array{Int64,2} = sum(A, dims=2)
    s::Float64 = sum(d)

    # Initialization
    f = vec(1 .+ 2 * d - (d .^ 2) / s)

    i::Int64 = 1
    diff::Float64 = 10.0

    Y0::Array{Float64,2} = zeros(size(X))
    Y::Array{Float64,2} = zeros(size(X))

    obj_0::Float64 = 1.0
    obj::Float64 = 1.0
    r::Float64 = 0.0

    while (i <= n_it) && (diff > t || i < 4)

        # Matrix product
        Y = A * X - (d / s) * (d' * X) + Diagonal(f) * X
        # In order to avoid calculating twice the same products
        # Since X is feasible, we calculate the objective value AFTER the first matrix product.
        obj = (tr(X' * Y)) / s

        #Difference between consecutive objectives
        diff = abs.(obj - obj_0) / obj_0

        #Update with momentum
        r = (i - 1.0) / (i + 2.0)
        X = (1 + r) * Y
        X = X - r * Y0

        #Project on the sphere
        X = X./sqrt.(sum(abs2, X, dims=2))

        #Update objective and momentum
        obj_0 = obj
        Y0 = Y
        i += 1
    end

    # Output
    if diff < t
        println("The iteration has become stationary after ", i,  " iterations")
    else
        println("The iteration did not converge after ",i," iterations")
        println("The relative difference between the last objective values", diff)
    end


    return X
end
