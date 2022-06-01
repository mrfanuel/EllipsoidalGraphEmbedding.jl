
function find_community_membership(x_embed::AbstractArray{Float64,2}, R0::Array{Float64,2})

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



function partition(A::SparseMatrixCSC{Int64,Int64}, x_embed::AbstractArray{Float64,2}, it_max::Int64, n_clu::Int64, p::Array{Float64,1})

    # Initialization
    d::Array{Int64,2} = sum(A, dims=2)
    s::Int64 = sum(d)

    #Number of nodes
    N::Int64 = size(x_embed, 2)

    # Random points as seeds for communities
    #index = random_generate_vec(p, n_clu) #index = rand(1:N,n_clu);
    index = rand(Categorical(p), n_clu)

    R0 = x_embed[:, index]
    dim::Int64 = size(x_embed, 1)

    ## Construct the initial communities and centroid vectors

    # Here is the first iteration
    community = find_community_membership(x_embed, R0)# Gives as many communities as nodes

    # Then, we update the centroids
    R = zeros(Float64, dim, n_clu)
    R, community = update_centroids(x_embed, R0, community)


    ## Update until stationarity

    community1 = vec(zeros(Int64, N, 1))
    R1 = zeros(dim, n_clu)

    n_updates::Int64 = 0
    community = rename_com_unique(vec(community))
    n_c = length(unique(community))

    H_lab = sparse(1:N, community, vec(ones(Int64, N, 1)), N, n_c)
    Q_best = (1 / s) * (tr(H_lab' * A * H_lab) - (norm(d' * H_lab, 2)^2) / s)

    for _ = 1:it_max
        R1, community1 = update_centroids(x_embed, R, community)
        community1 = rename_com_unique(vec(community1))
        n_c = length(unique(community1))
        H_lab = sparse(1:N, community1, vec(ones(Int64, N, 1)), N, n_c)
        Q = (1 / s) * (tr(H_lab' * A * H_lab) - (norm(d' * H_lab, 2)^2) / s)

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

function sphere_embed_cluster(A::SparseMatrixCSC{Int64,Int64}, n_it_PPM::Int64, t::Float64, n_clusters::Int64, n_it_Kmeans::Int64, n_updates::Int64, shape::String, r0::Int64)

    N = size(A, 1)
    d = sum(A, dims=2)
    s = sum(d)

    ############################# ROW NORMALIZATION ###################################################
    index::Array{Int64,1} = sample(1:N, r0; replace=false, ordered=true)
    H0 = A[:, index] - (1 / s) * d * d[index]'

    #Project on the sphere
    H0 = H0./sqrt.(sum(abs2, H0, dims=2));


    println(" ------- Acc Projected Power Iteration -------")
    fast::Int64 = 1
    H0 = @time acc_proj_power_method(A, H0, n_it_PPM, t, fast)
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

    for _ = 1:n_it_Kmeans
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



function update_centroids(S::AbstractArray{Float64,2}, R0::Array{Float64,2}, community::Array{Int64,1})

    N::Int64 = size(S, 2)
    dim::Int64 = size(S, 1)
    n_clu::Int64 = size(R0, 2)
    community1 = zeros(Int64, N, 1)

    # Update the centroids
    R1 = zeros(dim, n_clu)
    id::Int64 = 0

    for i = 1:N

        Si = view(S, :, i)
        si = R0' * Si
        si[community[i]] -= sum(Si .^ 2) # Correct for the "within the group effect"

        _, id = findmax(si)
        @inbounds community1[i] = id
        @inbounds R1id = view(R1, :, id)

        for j = 1:dim
            @inbounds R1id[j] += Si[j]
        end
    end

    return R1, community1
end

# function Weighted_Adjacency_List(w_edge_list::Array{Int64,2})

#     node_names = unique(w_edge_list[:, 1])
#     N = length(node_names)

#     n_connected_links = zeros(Int64, N, 1)
#     degree = zeros(Int64, N, 1)

#     node_names = unique(w_edge_list[:, 1])
#     n_connections = size(w_edge_list, 1)
#     for l = 1:n_connections
#         i = w_edge_list[l, 1]
#         index_i = find(node_names .== i)[1]
#         w = w_edge_list[l, 3]

#         n_connected_links[index_i] += 1
#         degree[index_i] += w
#     end

#     cumul_connected_links = zeros(Int64, N, 1)

#     cumul_connected_links[1] = 1
#     for index_l = 2:N
#         cumul_connected_links[index_l] = cumul_connected_links[index_l-1] + n_connected_links[index_l-1]
#     end

#     # Warning: The fact that the nodes are named from 1 to N is important
#     if sum(n_connected_links) !== n_connections
#         println("error")
#     end

#     w_adjacency_list = zeros(Int64, n_connections, 3)
#     for l = 1:n_connections
#         i = w_edge_list[l, 1]
#         index_i = find(node_names .== i)[1]
#         j = w_edge_list[l, 2]
#         w = w_edge_list[l, 3]

#         i_start = cumul_connected_links[index_i]
#         i_end = i_start + n_connected_links[index_i] - 1

#         ok = 0
#         k = i_start
#         while k <= i_end && ok == 0
#             if w_adjacency_list[k] == 0
#                 w_adjacency_list[k, 1] = i
#                 w_adjacency_list[k, 2] = j
#                 w_adjacency_list[k, 3] = w
#                 ok = 1
#             end
#             k += 1
#         end

#     end

#     return w_adjacency_list, n_connected_links, cumul_connected_links, degree

# end



# function ConComponentsCommunities(A::SparseMatrixCSC{Int64,Int64}, community::Array{Int64,1})

#     N = size(community, 1)
#     community_n = vec(zeros(Int64, N, 1))
#     n_c = length(unique(community))
#     # Find the largest component within each cluster
#     for i = 1:n_c
#         index_id = find(community .== i)
#         _, p1 = largest_component(A[index_id, index_id])
#         community_n[index_id[p1]] = i

#     end
#     return community_n
# end


# function FuseGraphEdgeList(w_edge_list::Array{Int64,2}, community::Array{Int64,1})

#     w_adjacency_list, n_connected_links, cumul_connected_links, degree = Weighted_Adjacency_List(w_edge_list::Array{Int64,2})


#     s = sum(degree)
#     N = maximum(size(degree))
#     n_connections = size(w_adjacency_list, 1)

#     #Construct representative nodes
#     community_names = unique(community)
#     N_com = countnz(community_names)
#     index_node_rep = zeros(Int64, N_com)

#     for i = 1:N_com
#         found = 0
#         j = 1
#         while found == 0
#             if community[j] == community_names[i]
#                 index_node_rep[i] = j
#                 found = 1
#             else
#                 j += 1
#             end
#         end
#     end


#     #Construct New EdgeList
#     edge_list_n = zeros(Int64, n_connections, 3)
#     l = 1
#     for index_i = 1:N
#         i_start = cumul_connected_links[index_i]
#         i_end = i_start + n_connected_links[index_i] - 1

#         # check if index_i is a representative 
#         if isempty(find(index_node_rep .== index_i))
#             index_i_n = index_node_rep[community[index_i]]

#         else
#             index_i_n = index_i
#         end

#         for k = i_start:i_end
#             index_j_k = w_adjacency_list[k, 2]

#             if isempty(find(index_node_rep .== index_j_k))
#                 index_j_n = index_node_rep[community[index_j_k]]
#             else
#                 index_j_n = index_j_k
#             end
#             wij = w_adjacency_list[k, 3]
#             edge_list_n[l, :] = [index_i_n index_j_n wij]

#             l += 1
#         end
#     end

#     # remove multiple lines
#     edge_list = unique(edge_list_n[:, 1:2], 1)
#     edge_list = [edge_list zeros(Int64, size(edge_list, 1), 1)]

#     for l = 1:size(edge_list, 1)

#         for k = 1:n_connections
#             if edge_list_n[k, 1] == edge_list[l, 1] && edge_list_n[k, 2] == edge_list[l, 2]
#                 edge_list[l, 3] += edge_list_n[k, 3]
#             end
#         end

#     end

#     community_n = zeros(Int64, N_com)

#     community_n = community[index_node_rep]

#     return edge_list, community_n

# end
# function GetElistFromAdjMatrix(A::SparseMatrixCSC{Int64,Int64})

#     #Construct a weighted edgelist from A
#     c = findnz(A)
#     n_connections = length(c[3])
#     w_edge_list = zeros(Int64, n_connections, 3)
#     for i = 1:n_connections
#         w_edge_list[i, 1] = c[2][i]
#         w_edge_list[i, 2] = c[1][i]
#         w_edge_list[i, 3] = c[3][i]
#     end
#     c = 0
#     gc()

#     return w_edge_list

# end


# function ConnectedCommunities(A::SparseMatrixCSC{Int64,Int64}, community::Array{Int64,1})
#     N = size(community, 1)
#     community_n = vec(zeros(Int64, N, 1))
#     n_c = length(unique(community))
#     # Find the largest component within each cluster
#     for i = 1:n_c
#         index_id = find(community .== i)
#         _, p1 = largest_component(A[index_id, index_id])
#         community_n[index_id[p1]] = i

#     end
#     # Put all the non-connected nodes in their own community
#     l = n_c + 1
#     for i = 1:N
#         if community_n[i] == 0
#             community_n[i] = l
#             l += 1
#         end
#     end

#     return community_n

# end


# function Convert2AdjMatrix(N::Int64, adj)
#     dim = size(adj)
#     N = maximum([maximum(adj[:, 1]), maximum(adj[:, 2])])
#     A = spzeros(N, N)
#     A = sparse(adj[:, 1], adj[:, 2], vec(ones(Int64, dim[1], 1)), N, N)
#     A0 = A + A'
#     A, _ = largest_component(A0)
#     return A
# end