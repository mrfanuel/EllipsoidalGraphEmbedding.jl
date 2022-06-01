module SphericalGraphEmbedding

using Distributions
using SparseArrays
using LinearAlgebra
using MatrixNetworks
using BenchmarkTools

include("utils.jl")
include("graph.jl")
include("optim.jl")


export
    random_generate,
    random_generate_vec,
    rename_com_unique,
    partition,
    sphere_embed_cluster,
    update_centroids,
    acc_proj_power_method
end
