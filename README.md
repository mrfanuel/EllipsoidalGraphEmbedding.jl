# EllipsoidalGraphEmbedding.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mrfanuel.github.io/SphericalGraphEmbedding.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mrfanuel.github.io/SphericalGraphEmbedding.jl/dev)
[![Build Status](https://github.com/mrfanuel/SphericalGraphEmbedding.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mrfanuel/SphericalGraphEmbedding.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mrfanuel/SphericalGraphEmbedding.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mrfanuel/SphericalGraphEmbedding.jl)


### Documentation
This is the code associated with the manuscript 
[Ellipsoidal Embedding of graphs](https://github.com/mrfanuel/SphericalGraphEmbedding.jl)
by [Michaël Fanuel](https://mrfanuel.github.io/) , [Antoine Aspeel](https://scholar.google.com/citations?user=EDDQMfgAAAAJ&hl=en), [Michael Schaub](https://michaelschaub.github.io/) and [Jean-Charles Delvenne](https://perso.uclouvain.be/jean-charles.delvenne/welcome.html).

### Install Julia

If you do not have Julia installed, please visit [julialang.org](https://julialang.org/learning/getting-started/)
### Installation

[`EllipsoidalGraphEmbedding.jl`](https://github.com/mrfanuel/EllipsoidalGraphEmbedding.jl) is not registered.
The way to use it is to type

```julia
julia> ]add https://github.com/mrfanuel/EllipsoidalGraphEmbedding.jl
```

### Jupyter notebooks for reproducing the paper figures

You can execute the Jupyter [`notebooks`](https://github.com/mrfanuel/EllipsoidalGraphEmbedding.jl/blob/master/notebooks) to generate the paper figures.

- [spherical_embedding_plots.ipynb](https://github.com/mrfanuel/EllipsoidalGraphEmbedding.jl/blob/main/notebooks/spherical_embedding_plots.ipynb) reproduces Figure 1.1, Figure 1.2 and Figure 4.1.
- [nmi_analysis.ipynb](https://github.com/mrfanuel/EllipsoidalGraphEmbedding.jl/blob/main/notebooks/nmi_analysis.ipynb) yields Figure 4.2 and Figure C.1.
- [ppm.ipynb](https://github.com/mrfanuel/EllipsoidalGraphEmbedding.jl/blob/main/notebooks/ppm.ipynb) gives Figure B.1.
- [real_networks_benchmark.ipynb](https://github.com/mrfanuel/EllipsoidalGraphEmbedding.jl/blob/main/notebooks/real_networks_benchmark.ipynb) produces results of ellipsoidal and spectral embeddings in Table C.2.

### Usage

To use this functions of this package, simply type

```julia
julia> using EllipsoidalGraphEmbedding
```
