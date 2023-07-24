"""
    find_community_membership(x_embed, R0)

returns the index of the closet centroid to each position vector.
# Arguments
- `x_embed::AbstractArray{Float64,2}` array of position vectors
- `R0::Array{Float64,2}`array of position vectors for centroids.
# Output
- `community::Array{Int64,1}` membership array

"""
function find_community_membership(x_embed::AbstractArray{Float64,2}, R0::Array{Float64,2})#::Array{Int64,1}

    N::Int64 = size(x_embed, 2)
    id::Int64 = 0
    community = vec(zeros(Int64, 1, N))

    for i = 1:N
        @inbounds si = R0' * x_embed[:, i]
        _, id = findmax(si)
        @inbounds community[i] = id
    end

    return community

end


@doc raw"""
partition(A, x_embed, it_max, n_clusters, p)

returns a clustering of the embedded nodes.
# Arguments
- `A::SparseMatrixCSC{Int64,Int64}` adjacency matrix
- `x_embed::AbstractArray{Float64,2}` array of position vectors
- `it_max::Int64` maximum number of iterations
- `n_clusters::Int64` maximum number of clusters
- `p::Array{Float64,1}` array of probabilities for sampling the centroids

# Output
- `community::Array{Int64,1}` membership array
- `Q_best::Float64` modularity of the partition
- `n_updates::Float64` number of centroid updates to reach stationarity

"""
function partition(A::SparseMatrixCSC{Int64,Int64}, x_embed::AbstractArray{Float64,2}, it_max::Int64, n_clusters::Int64, p::Array{Float64,1})#::Tuple{Array{Int64,1},Float64,Int64}

    # Initialization
    d::Array{Int64,2} = sum(A, dims=2)
    s::Int64 = sum(d)

    #Number of nodes
    N::Int64 = size(x_embed, 2)

    # Random points as seeds for communities
    index = rand(Categorical(p), n_clusters)

    R0 = x_embed[:, index]
    dim::Int64 = size(x_embed, 1)

    ## Construct the initial communities and centroid vectors

    # Here is the first iteration
    community = find_community_membership(x_embed, R0)# Gives as many communities as nodes

    # Then, we update the centroids
    R = zeros(Float64, dim, n_clusters)
    R, community = update_centroids(x_embed, R0, community)


    ## Update until stationarity

    community1 = vec(zeros(Int64, N, 1))
    R1 = zeros(dim, n_clusters)

    n_updates::Int64 = 0
    community = rename_com_unique(vec(community))
    n_c = length(unique(community))

    H_lab = sparse(1:N, community, vec(ones(Int64, N, 1)), N, n_c)

    Q_best = tr(H_lab' * x_embed' * x_embed * H_lab)

    for _ = 1:it_max
        R1, community1 = update_centroids(x_embed, R, community)
        community1 = rename_com_unique(vec(community1))
        n_c = length(unique(community1))
        H_lab = sparse(1:N, community1, vec(ones(Int64, N, 1)), N, n_c)

        Q = tr(H_lab' * x_embed' * x_embed * H_lab)

        if Q > Q_best
            R = R1
            community = community1
            Q_best = Q
            n_updates = n_updates + 1
        else
            break
        end
    end

    return community, Q_best, n_updates
end

@doc raw"""
    sphere_embed_cluster(A, n_it_PPM, t, n_clusters, n_rep_vec_part, n_updates, shape, r0)
returns a clustering of the embedded nodes.
# Arguments
- `A::SparseMatrixCSC{Int64,Int64}` adjacency matrix
- `n_it_PPM::Int64` maximum number of iterations for project power method (PPM)
- `tol::Float64` tolerance for relative objective variation for PPM
- `n_clusters::Int64` maximum number of clusters
- `n_rep_vec_part::Int64` number of repetitions of vector partition
- `n_updates::Int64` number of updates of vector partition
- `shape::String` spherical or ellipsoidal embedding
- `r0::Int64` largest number of dimensions for the embedding

# Output
- `x_embed::AbstractArray{Float64,2}` array of position vectors
- `community::Array{Int64,1}` membership array

"""
function sphere_embed_cluster(A::SparseMatrixCSC{Int64,Int64}, n_it_PPM::Int64, tol::Float64, n_clusters::Int64, n_rep_vec_part::Int64, n_updates::Int64, shape::String, r0::Int64)#::Tuple{AbstractArray{Float64,2},Array{Int64,1}}

    N = size(A, 1)
    d = sum(A, dims=2)
    s = sum(d)

    ############################# ROW NORMALIZATION ###################################################
    index::Array{Int64,1} = sample(1:N, r0; replace=false, ordered=true)
    H0 = A[:, index] - (1 / s) * d * d[index]'

    #Project on the sphere
    H0 = H0./sqrt.(sum(abs2, H0, dims=2));


    println(" ------- Acc Projected Power Iteration -------")
    H0 = @time acc_proj_power_method(A, H0, n_it_PPM, tol)
    U, S, _ = svd(H0)
    println(" -------------- Clustering ------- ")
    dim::Int64 = length(S)
    println("dimension of embedding used for clustering: ",dim)

    x_embed = zeros(dim, N)
    if shape == "Spherical"
        x_embed = (U[:, 1:dim] * diagm((S[1:dim])))'
    elseif shape == "Ellipsoidal"
        x_embed = (U[:, 1:dim])'
    end

    Q_best::Float64 = 0
    n_c_best::Int64 = 0
    n_c::Int64 = 0

    # Probability mass for sampling the initial centroids
    p = vec(d / s)
    n_updates_best::Int64 = 0

    ## keep partition with best modularity

    # initialization
    community, Q, n_updates = partition(A, x_embed, n_updates, n_clusters, p)

    for _ = 1:n_rep_vec_part
        community0, Q, n_updates = partition(A, x_embed, n_updates, n_clusters, p)
        n_c = length(unique(community0))
        if Q > Q_best
            community = community0
            Q_best = Q
            n_c_best = n_c
            n_updates_best = n_updates
        end
    end

    community = rename_com_unique(community)

    print("Number of updates: ")
    println(n_updates_best)

    print("Number of communities: ")
    println(n_c_best)

    print("Modularity: ")
    println(Q_best)
    println(" -------------------------------------------- ")
    println("The first 5 squared singular values divided by N : ")
    println((S[1:5] .^ 2) / N)
    println(" -------------------------------------------- ")

    return x_embed, community

end

"""
    update_centroids(S, R0, community)

returns a clustering of the embedded nodes.
# Arguments
- `S::AbstractArray{Float64,2}` position vectors of embedded nodes
- `R0::Array{Float64,2}` initial position of the centroids
- `community::Array{Int64,1}` initial node membership array 


# Output
- `R1::Array{Float64,2}` updated position of the centroids
- `community1::Array{Int64,1}` updated node membership array 

"""
function update_centroids(S::AbstractArray{Float64,2}, R0::Array{Float64,2}, community0::Array{Int64,1})

    N::Int64 = size(S, 2)
    dim::Int64 = size(S, 1)
    n_clusters::Int64 = size(R0, 2)
    community1 = zeros(Int64, N, 1)

    # Update the centroids
    R1 = zeros(dim, n_clusters)
    id::Int64 = 0

    for i = 1:N

        Si = view(S, :, i)
        si = R0' * Si
        si[community0[i]] -= sum(Si .^ 2) # Correct for the "within the group effect"

        _, id = findmax(si)
        @inbounds community1[i] = id
        @inbounds R1id = view(R1, :, id)

        for j = 1:dim
            @inbounds R1id[j] += Si[j]
        end
    end

    return R1, community1
end

