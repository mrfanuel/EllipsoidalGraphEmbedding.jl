var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = SphericalGraphEmbedding","category":"page"},{"location":"#SphericalGraphEmbedding","page":"Home","title":"SphericalGraphEmbedding","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for SphericalGraphEmbedding.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [SphericalGraphEmbedding]","category":"page"},{"location":"#SphericalGraphEmbedding.acc_proj_power_method-Tuple{SparseArrays.SparseMatrixCSC{Int64, Int64}, Matrix{Float64}, Int64, Float64}","page":"Home","title":"SphericalGraphEmbedding.acc_proj_power_method","text":"acc_proj_power_method(A, X, n_it, tol)\n\nRuns the accelarated projected power iteration on the modularity matrix.\n\nArguments\n\nA::SparseMatrixCSC{Int64,Int64} sparse graph adjacency matrix.\nX::Array{Float64,2} initial matrix.\nn_it::Int64 maximum number of iterations.\ntol::Float64 tolerance on relative variation of consecutive objective values.\n\nOutput\n\nX::Array{Float64,2} approximate stationary matrix of the iteration.\n\n\n\n\n\n","category":"method"},{"location":"#SphericalGraphEmbedding.find_community_membership-Tuple{AbstractMatrix{Float64}, Matrix{Float64}}","page":"Home","title":"SphericalGraphEmbedding.find_community_membership","text":"find_community_membership(x_embed, R0)\n\nreturns the index of the closet centroid to each position vector.\n\nArguments\n\nx_embed::AbstractArray{Float64,2} array of position vectors\nR0::Array{Float64,2}array of position vectors for centroids.\n\nOutput\n\ncommunity::Array{Int64,1} membership array\n\n\n\n\n\n","category":"method"},{"location":"#SphericalGraphEmbedding.partition-Tuple{SparseArrays.SparseMatrixCSC{Int64, Int64}, AbstractMatrix{Float64}, Int64, Int64, Vector{Float64}}","page":"Home","title":"SphericalGraphEmbedding.partition","text":"partition(A, xembed, itmax, n_clusters, p)\n\nreturns a clustering of the embedded nodes.\n\nArguments\n\nA::SparseMatrixCSC{Int64,Int64} adjacency matrix\nx_embed::AbstractArray{Float64,2} array of position vectors\nit_max::Int64 maximum number of iterations\nn_clusters::Int64 maximum number of clusters\np::Array{Float64,1} array of probabilities for sampling the centroids\n\nOutput\n\ncommunity::Array{Int64,1} membership array\nQ_best::Float64 modularity of the partition\nn_updates::Float64 number of centroid updates to reach stationarity\n\n\n\n\n\n","category":"method"},{"location":"#SphericalGraphEmbedding.rename_com_unique-Tuple{Vector{Int64}}","page":"Home","title":"SphericalGraphEmbedding.rename_com_unique","text":"rename_com_unique(community)\n\nrename the community labels by the smallest labels of integers.\n\nArguments\n\ncommunity::Array{Int64,1} array of integers\n\nOutput\n\ncommunity::Array{Int64,1} renamed array of integers\n\n\n\n\n\n","category":"method"},{"location":"#SphericalGraphEmbedding.spectral_embed_cluster-Tuple{SparseArrays.SparseMatrixCSC{Int64, Int64}, Int64, Float64, Vararg{Int64, 4}}","page":"Home","title":"SphericalGraphEmbedding.spectral_embed_cluster","text":"sphere_embed_cluster(A, n_it_PPM, t, n_clusters, n_rep_vec_part, n_updates, shape, r0)\n\nreturns a clustering of the embedded nodes.\n\nArguments\n\nA::SparseMatrixCSC{Int64,Int64} adjacency matrix\nit_max::Int64 maximum number of iterations of power method with Gram Schmidt\ntol::Float64 tolerance for power method\nn_clusters::Int64 maximum number of clusters\nn_rep_vec_part::Int64 number of repetitions of vector partition\nn_updates::Int64 number of updates of vector partition\ndim_embed_spectral::Int64 largest number of dimensions for the embedding\n\nOutput\n\nx_embed::AbstractArray{Float64,2} array of position vectors\ncommunity::Array{Int64,1} membership array\nS::Vector{Float64} eigenvalues spectral embedding\n\n\n\n\n\n","category":"method"},{"location":"#SphericalGraphEmbedding.sphere_embed_cluster-Tuple{SparseArrays.SparseMatrixCSC{Int64, Int64}, Int64, Float64, Int64, Int64, Int64, String, Int64}","page":"Home","title":"SphericalGraphEmbedding.sphere_embed_cluster","text":"sphere_embed_cluster(A, n_it_PPM, t, n_clusters, n_rep_vec_part, n_updates, shape, r0)\n\nreturns a clustering of the embedded nodes.\n\nArguments\n\nA::SparseMatrixCSC{Int64,Int64} adjacency matrix\nn_it_PPM::Int64 maximum number of iterations for project power method (PPM)\ntol::Float64 tolerance for relative objective variation for PPM\nn_clusters::Int64 maximum number of clusters\nn_rep_vec_part::Int64 number of repetitions of vector partition\nn_updates::Int64 number of updates of vector partition\nshape::String spherical or ellipsoidal embedding\nr0::Int64 largest number of dimensions for the embedding\n\nOutput\n\nx_embed::AbstractArray{Float64,2} array of position vectors\ncommunity::Array{Int64,1} membership array\nS::Vector{Float64} singular values embedding\n\n\n\n\n\n","category":"method"},{"location":"#SphericalGraphEmbedding.update_centroids-Tuple{AbstractMatrix{Float64}, Matrix{Float64}, Vector{Int64}}","page":"Home","title":"SphericalGraphEmbedding.update_centroids","text":"update_centroids(S, R0, community)\n\nreturns a clustering of the embedded nodes.\n\nArguments\n\nS::AbstractArray{Float64,2} position vectors of embedded nodes\nR0::Array{Float64,2} initial position of the centroids\ncommunity::Array{Int64,1} initial node membership array \n\nOutput\n\nR1::Array{Float64,2} updated position of the centroids\ncommunity1::Array{Int64,1} updated node membership array \n\n\n\n\n\n","category":"method"}]
}
