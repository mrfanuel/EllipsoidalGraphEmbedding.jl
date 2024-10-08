module EllipsoidalGraphEmbedding

using Distributions
using SparseArrays
using LinearAlgebra
using MatrixNetworks
using BenchmarkTools
using Random

include("utils.jl")
include("graph.jl")
include("optim.jl")


export
    rename_com_unique,
    partition,
    sphere_embed_cluster,
    update_centroids,
    acc_proj_power_method,
    mat_vec_prod,
    gram_schmidt,
    top_eigenpairs_Q,
    spectral_embed_cluster,
    dim_eff,
    adj_matrix_cc,
    largest_cc_from_adj
end