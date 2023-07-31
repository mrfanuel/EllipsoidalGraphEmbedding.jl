"""
    rename_com_unique(community)

rename the community labels by the smallest labels of integers.
# Arguments
- `community::Array{Int64,1}` array of integers

# Output
- `community::Array{Int64,1}` renamed array of integers

"""
function rename_com_unique(community::Array{Int64,1})#::Array{Int64,1}

    N = length(community)
    community_reordered = vec(zeros(Int64, N, 1))

    n_c_best = length(unique(community))
    index = unique(community)
    for i = 1:N
        for l = 1:n_c_best
            if community[i] == index[l]
                community_reordered[i] = l
            end
        end
    end
    return community_reordered
end


function mat_vec_prod(A,v0,d,shift)
    sum_d = sum(d)
    v0 = A * v0 .- (d / sum_d) * (d' * v0) .+ shift * v0
    return v0
end

# since I could not find a thin qr decomposition, simply use Gram Schmidt process
function gram_schmidt(A)
	# performs Gram Schmidt process on columns of matrix A
	# the columns of the returned matrix Q are the outputs of this process
	# https://www.math.hkust.edu.hk/~emarberg/teaching/2021/Math2121/julia/14_Math2121_Fall2021.jl.html
    m = size(A, 1)
	n = size(A, 2)
    Q = zeros(m, n)
    norm_col = zeros(n,1)
    
    for j=1:n
        v = A[1:m, j]
		x = A[1:m, j]
        for i = 1:j - 1
            u = Q[1:m, i]
            x -= dot(u, v) * u / dot(u, u)
        end
        t = norm(x)
        norm_col[j] = t
        Q[1:m, j] =  x / t
    end

    return Q , norm_col
end

function top_eigenpairs_Q(A,r,tol,it_max)
    # compute r leading eigenpairs of modularity matrix Q = A - dd'/s
    # by avoiding construction of Q

    N = size(A)[1]
    d = sum(A,dims=2)
    sum_d = sum(d)

    # Preprocessing: compute shift, i.e.,
    # an estimate of absolute largest eigenvalue
    # to make Q + shift * I positive semidefinite

    v0 = rand(Uniform(-1,1),N,1)
    v0 = v0 / norm(v0)

    n=0
    shift = 0
    for it = 1:30
        v0 = mat_vec_prod(A,v0,d,shift)
        n = norm(v0)
        v0 = v0 / n
    end
    shift = n

    # initialization
    v0 = rand(Uniform(-1,1),N,r)
    v0 = v0 ./ sqrt.(sum(abs2, v0, dims=1))
    v0 , n_col0 = gram_schmidt(v0)

    # Krylov method for r top eigenvalues
    it = 1
    while it <= it_max
        v0 = mat_vec_prod(A,v0,d,shift)
        v0 , n_col = gram_schmidt(v0)
        it += 1
        if norm(n_col - n_col0) / norm(n_col0) < tol
            break
        end
        n_col0 = n_col
    end

    if it == it_max
        println("not converged")
    end

    eigenvalues = n_col0 .- shift

    return eigenvalues , v0
end

function dim_eff(singular,epsilon)

    s = 0
    s_tot = sum(singular.^2)
    d_eff = length(singular)
    for i = 1:length(singular)
        s += singular[i]^2
        p = s/s_tot
        if p > 1 - epsilon 
            d_eff = i
            break
        end
    end
    return d_eff
end