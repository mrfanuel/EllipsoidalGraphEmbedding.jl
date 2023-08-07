var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = EllipsoidalGraphEmbedding","category":"page"},{"location":"#EllipsoidalGraphEmbedding","page":"Home","title":"EllipsoidalGraphEmbedding","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for EllipsoidalGraphEmbedding.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [EllipsoidalGraphEmbedding]","category":"page"},{"location":"#EllipsoidalGraphEmbedding.tol","page":"Home","title":"EllipsoidalGraphEmbedding.tol","text":"sphere_embed_cluster(A, n_it_PPM, t, n_clusters, n_rep_vec_part, n_updates, shape, r0)\n\nreturns a clustering of the embedded nodes.\n\nArguments\n\nA::SparseMatrixCSC{Int64,Int64} adjacency matrix\ndim_embed_spectral::Int64 largest number of dimensions for the embedding\nn_clusters::Int64 maximum number of clusters\nit_max::Int64 maximum number of iterations of power method with Gram Schmidt\ntol::Float64 tolerance for power method\nn_rep_vec_part::Int64 number of repetitions of vector partition\nn_updates::Int64 number of updates of vector partition\n\nOutput\n\nx_embed::AbstractArray{Float64,2} array of position vectors\ncommunity::Array{Int64,1} membership array\nS::Vector{Float64} eigenvalues spectral embedding\n\n\n\n\n\n","category":"constant"},{"location":"#EllipsoidalGraphEmbedding.acc_proj_power_method","page":"Home","title":"EllipsoidalGraphEmbedding.acc_proj_power_method","text":"acc_proj_power_method(A, X, n_it, tol, descriptor)\n\nRuns the accelarated projected power iteration on the modularity matrix.\n\nArguments\n\nA::SparseMatrixCSC{Int64,Int64} sparse graph adjacency matrix.\nX::Array{Float64,2} initial matrix.\nn_it::Int64 maximum number of iterations.\ntol::Float64 tolerance on relative variation of consecutive objective values.\n`descriptor::String' type of descriptor matrix, by default \"Modularity\", otherwise \"Laplacian\".\n\nOutput\n\nX::Array{Float64,2} approximate stationary matrix of the iteration.\n\n\n\n\n\n","category":"function"},{"location":"#EllipsoidalGraphEmbedding.find_community_membership-Tuple{AbstractMatrix{Float64}, Matrix{Float64}}","page":"Home","title":"EllipsoidalGraphEmbedding.find_community_membership","text":"find_community_membership(x_embed, R0)\n\nreturns the index of the closet centroid to each position vector.\n\nArguments\n\nx_embed::AbstractArray{Float64,2} array of position vectors\nR0::Array{Float64,2}array of position vectors for centroids.\n\nOutput\n\ncommunity::Array{Int64,1} membership array\n\n\n\n\n\n","category":"method"},{"location":"#EllipsoidalGraphEmbedding.partition-Tuple{AbstractMatrix{Float64}, Int64, Int64, Vector{Float64}}","page":"Home","title":"EllipsoidalGraphEmbedding.partition","text":"partition(xembed, itmax, n_clusters, p)\n\nreturns a clustering of the embedded nodes.\n\nArguments\n\nx_embed::AbstractArray{Float64,2} array of position vectors\nit_max::Int64 maximum number of iterations\nn_clusters::Int64 maximum number of clusters\np::Array{Float64,1} array of probabilities for sampling the centroids\n\nOutput\n\ncommunity::Array{Int64,1} membership array\nObj_best::Float64 objective of the best partition\nn_updates::Float64 number of centroid updates to reach stationarity\n\n\n\n\n\n","category":"method"},{"location":"#EllipsoidalGraphEmbedding.rename_com_unique-Tuple{Vector{Int64}}","page":"Home","title":"EllipsoidalGraphEmbedding.rename_com_unique","text":"rename_com_unique(community)\n\nrename the community labels by the smallest labels of integers.\n\nArguments\n\ncommunity::Array{Int64,1} array of integers\n\nOutput\n\ncommunity::Array{Int64,1} renamed array of integers\n\n\n\n\n\n","category":"method"},{"location":"#EllipsoidalGraphEmbedding.sphere_embed_cluster","page":"Home","title":"EllipsoidalGraphEmbedding.sphere_embed_cluster","text":"sphere_embed_cluster(A, d0, shape,n_clusters, n_it_PPM, tol,  n_rep_vec_part, n_updates,descriptor)\n\nreturns a clustering of the embedded nodes.\n\nArguments\n\nA::SparseMatrixCSC{Int64,Int64} adjacency matrix\nd0::Int64 largest number of dimensions for the embedding\nshape::String \"Spherical\" or \"Ellipsoidal\"\nn_clusters::Int64 maximum number of clusters\nn_it_PPM::Int64 maximum number of iterations for project power method (PPM)\ntol::Float64 tolerance for relative objective variation for PPM\nn_rep_vec_part::Int64 number of repetitions of vector partition\nn_updates::Int64 number of updates of vector partition\ndescriptor::String = \"Modularity\" or \"Laplacian\"\n\nOutput\n\nx_embed::AbstractArray{Float64,2} array of position vectors\ncommunity::Array{Int64,1} membership array\nS::Vector{Float64} singular values embedding\n\n\n\n\n\n","category":"function"},{"location":"#EllipsoidalGraphEmbedding.update_centroids-Tuple{AbstractMatrix{Float64}, Matrix{Float64}, Vector{Int64}}","page":"Home","title":"EllipsoidalGraphEmbedding.update_centroids","text":"update_centroids(S, R0, community)\n\nreturns a clustering of the embedded nodes.\n\nArguments\n\nS::AbstractArray{Float64,2} position vectors of embedded nodes\nR0::Array{Float64,2} initial position of the centroids\ncommunity::Array{Int64,1} initial node membership array \n\nOutput\n\nR1::Array{Float64,2} updated position of the centroids\ncommunity1::Array{Int64,1} updated node membership array \n\n\n\n\n\n","category":"method"}]
}
