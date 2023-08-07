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
partition(x_embed, it_max, n_clusters, p)

returns a clustering of the embedded nodes.
# Arguments
- `x_embed::AbstractArray{Float64,2}` array of position vectors
- `it_max::Int64` maximum number of iterations
- `n_clusters::Int64` maximum number of clusters
- `p::Array{Float64,1}` array of probabilities for sampling the centroids

# Output
- `community::Array{Int64,1}` membership array
- `Obj_best::Float64` objective of the best partition
- `n_updates::Float64` number of centroid updates to reach stationarity

"""
function partition(x_embed::AbstractArray{Float64,2}, it_max::Int64, n_clusters::Int64, p::Array{Float64,1})#::Tuple{Array{Int64,1},Float64,Int64}


    # Number of nodes
    N::Int64 = size(x_embed, 2)

    # Random points as seeds for communities
    index = rand(Categorical(p), n_clusters)

    R0 = x_embed[:, index]
    dim::Int64 = size(x_embed, 1)

    # Construct the initial communities and centroid vectors

    # Here is the first iteration
    community = find_community_membership(x_embed, R0)# Gives as many communities as nodes

    # Then, we update the centroids
    R = zeros(Float64, dim, n_clusters)
    R, community = update_centroids(x_embed, R0, community)


    # Update until stationarity

    community1 = vec(zeros(Int64, N, 1))
    R1 = zeros(dim, n_clusters)

    n_updates::Int64 = 0
    community = rename_com_unique(vec(community))
    n_c = length(unique(community))

    H_lab = sparse(1:N, community, vec(ones(Int64, N, 1)), N, n_c)

    Obj_best = tr(H_lab' * x_embed' * x_embed * H_lab)

    for _ = 1:it_max
        R1, community1 = update_centroids(x_embed, R, community)
        community1 = rename_com_unique(vec(community1))
        n_c = length(unique(community1))
        H_lab = sparse(1:N, community1, vec(ones(Int64, N, 1)), N, n_c)

        Obj = tr(H_lab' * x_embed' * x_embed * H_lab)

        if Obj > Obj_best
            R = R1
            community = community1
            Obj_best = Obj
            n_updates = n_updates + 1
        else
            break
        end
    end

    return community, Obj_best, n_updates
end

@doc raw"""
    sphere_embed_cluster(A, d0, shape,n_clusters, n_it_PPM, tol,  n_rep_vec_part, n_updates,descriptor)
returns a clustering of the embedded nodes.
# Arguments
- `A::SparseMatrixCSC{Int64,Int64}` adjacency matrix
- `d0::Int64` largest number of dimensions for the embedding
- `shape::String` "Spherical" or "Ellipsoidal"
- `n_clusters::Int64` maximum number of clusters

- `n_it_PPM::Int64` maximum number of iterations for project power method (PPM)
- `tol::Float64` tolerance for relative objective variation for PPM
- `n_rep_vec_part::Int64` number of repetitions of vector partition
- `n_updates::Int64` number of updates of vector partition
- `descriptor::String` = "Modularity" or "Laplacian"
# Output
- `x_embed::AbstractArray{Float64,2}` array of position vectors
- `community::Array{Int64,1}` membership array
- `S::Vector{Float64}` singular values embedding

"""
function sphere_embed_cluster(A::SparseMatrixCSC{Int64,Int64}, r0::Int64, shape::String,n_clusters::Int64, n_it_PPM::Int64 = 30000, tol::Float64 = 1e-8,  n_rep_vec_part::Int64 = 5, n_updates::Int64 = 50,descriptor::String = "Modularity")#::Tuple{AbstractArray{Float64,2},Array{Int64,1}}

    N = size(A, 1)
    d = sum(A, dims=2)
    s = sum(d)

    # ROW NORMALIZATION 
    index::Array{Int64,1} = sample(1:N, r0; replace=false, ordered=true)
    H0 = A[:, index] - (1 / s) * d * d[index]'

    # Project on the sphere
    H0 = H0./sqrt.(sum(abs2, H0, dims=2));


    println(" ------- Acc Projected Power Iteration -------")
    H0 = @time acc_proj_power_method(A, H0, n_it_PPM, tol, descriptor)
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

    Obj_best::Float64 = 0
    n_c_best::Int64 = 0
    n_c::Int64 = 0

    # Probability mass for sampling the initial centroids
    p = vec(d / s)
    n_updates_best::Int64 = 0

    # keep partition with best tr(H_lab' * x_embed' * x_embed * H_lab)

    # initialization
    community, Obj, n_updates = partition(x_embed, n_updates, n_clusters, p)

    for _ = 1:n_rep_vec_part
        community0, Obj, n_updates = partition(x_embed, n_updates, n_clusters, p)
        n_c = length(unique(community0))
        if Obj > Obj_best
            community = community0
            Obj_best = Obj
            n_c_best = n_c
            n_updates_best = n_updates
        end
    end
    print("Number of updates: ")
    println(n_updates_best)

    community = rename_com_unique(community)
    n_com = length(community)
    
    H_lab = sparse(1:N, community, vec(ones(Int64, N, 1)), N, n_com)
    modularity = (1 / s) * (tr(H_lab' * A * H_lab) - (norm(d' * H_lab, 2)^2) / s)
    print("Modularity: ")
    println(modularity)

    print("Number of communities: ")
    println(n_c_best)


    println(" -------------------------------------------- ")
    if r0 >= 5
        println("The first 5 squared singular values divided by N : ")
        println((S[1:5] .^ 2) / N)
    else
        println("The squared singular values divided by N : ")
        println((S .^ 2) / N)
    end
    println(" -------------------------------------------- ")

    return x_embed, community, S

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

@doc raw"""
    sphere_embed_cluster(A, n_it_PPM, t, n_clusters, n_rep_vec_part, n_updates, shape, r0)
returns a clustering of the embedded nodes.
# Arguments
- `A::SparseMatrixCSC{Int64,Int64}` adjacency matrix
- `dim_embed_spectral::Int64` largest number of dimensions for the embedding
- `n_clusters::Int64` maximum number of clusters
- `it_max::Int64` maximum number of iterations of power method with Gram Schmidt
- `tol::Float64` tolerance for power method
- `n_rep_vec_part::Int64` number of repetitions of vector partition
- `n_updates::Int64` number of updates of vector partition


# Output
- `x_embed::AbstractArray{Float64,2}` array of position vectors
- `community::Array{Int64,1}` membership array
- `S::Vector{Float64}` eigenvalues spectral embedding

"""
tol = 1e-06
it_max = 1000
function spectral_embed_cluster(A::SparseMatrixCSC{Int64,Int64}, dim_embed_spectral::Int64, n_clusters::Int64, it_max::Int64=1000, tol::Float64=1e-06, n_rep_vec_part::Int64= 5, n_updates::Int64=50)#::Tuple{AbstractArray{Float64,2},Array{Int64,1}}

    N = size(A, 1)
    d = sum(A, dims=2)

    eigenvalues , v0 = @time top_eigenpairs_Q(A,dim_embed_spectral,tol,it_max)

    d = sum(A,dims=2)
    sum_d = sum(d)
    println("spectral embedding computed")

    # possible test (passed for PowerEU)
    # using Arpack
    # Q = A - (1 / sum_d) * d * d' 
    # λ, ϕ = eigs(Q, nev = 10, which=:LR)
    # norm(eigenvalues - λ) / norm(eigenvalues)

    # Probability mass for sampling the initial centroids
    sum_d = sum(d)

    p = vec(d / sum_d)
    # initialization
    modularities = zeros(dim_embed_spectral,1)
    communities = zeros(dim_embed_spectral,N)
    nb_communities = zeros(dim_embed_spectral,1)

    for n_ev = 1:dim_embed_spectral
        Obj_best = 0
        n_c_best = 0
        n_c = 0

        n_updates_best = 0
        if n_ev > 1
            x_embed = (v0[:,1:n_ev])'

            # keep partition with tr(H_lab' * x_embed' * x_embed * H_lab)

            # initialization
            community, Obj, n_updates = partition(x_embed, n_updates, n_clusters, p)

            for _ = 1:n_rep_vec_part
                community0, Obj, n_updates = partition(x_embed, n_updates, n_clusters, p)
                n_c = length(unique(community0))
                if Obj > Obj_best
                    community = community0
                    Obj_best = Obj
                    n_c_best = n_c
                    n_updates_best = n_updates
                end
            end

            community = rename_com_unique(community)

            n_com = length(community)
            H_lab = sparse(1:N, community, vec(ones(Int64, N, 1)), N, n_com)
            modularity = (1 / sum_d) * (tr(H_lab' * A * H_lab) - (norm(d' * H_lab, 2)^2) / sum_d)

            modularities[n_ev] = modularity
            nb_communities[n_ev] = n_c
            communities[n_ev,:] = community
        end
    end

    mod_best, n_ev_best = findmax(modularities)

    x_embed = v0'
    community = communities[n_ev_best,:]
    r = vec(1:N)
    n_ev_best = r[n_ev_best]
    println("Number of eigenvector for best modularity: ", n_ev_best)
    println("Number of communities: ", nb_communities[n_ev_best])
    println("Modularity: ", mod_best)

    return x_embed, community, eigenvalues

end