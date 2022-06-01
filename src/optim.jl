# Input: 
#
# A : sparse adjacency matrix N x N
# H0 : feasible initial matrix N x r0
# n_it : maximum number of iterations
# t : tolerance (used for assessing the stationarity of the sequence of iterates).
# fast: if ==1 this option approximates the largest eigenvalue of the modularity matrix (for large matrices)
# fast: if ==0 this option calculates the largest eigenvalue of the modularity matrix (for small matrices)
#
# Output:
#
# X : a feasible matrix N x r0

function acc_proj_power_method(A::SparseMatrixCSC{Int64,Int64}, X::Array{Float64,2}, n_it::Int64, t::Float64, fast::Int64)

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
