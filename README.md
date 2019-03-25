# Topological Data Analysis

This package provides various tools for topological data analysis.

## Installation

```julia
pkg> add https://github.com/wildart/SmithNormalForm.jl.git#0.2.1
pkg> add https://github.com/wildart/ComputationalHomology.jl.git#master
pkg> add https://github.com/wildart/TDA.jl.git#master
```

For Julia 1.1+, add [BoffinStuff](https://github.com/wildart/BoffinStuff.git) registry in package manager, and proceed installation:

```
pkg> registry add https://github.com/wildart/BoffinStuff.git
pkg> add TDA
```

## Examples

### Persistance Diagram & Barcode

```julia
using TDA, ComputationalHomology, Plots
# crate some intervals of various dimensions
ints = vcat(intervals(0, 2.0=>6.0, 5.0=>10.0, 1.0=>Inf), intervals(1, 9.0=>12.0))
# plot persistance diagram
plot(ints)
# plot barcode
plot(ints, seriestype=:barcode)
```

### Nerve

```julia
using TDA, ComputationalHomology, Plots
# generate simplicial complex
cplx = ComputationalHomology.sphere(2)
# generate some points on circle
D = mapslices(p->p./sqrt(sum(p.^2)), randn(30,2), dims=2)
# plot points
plot(D[:,1], D[:,2], seriestype = :scatter, markersize = 2)
# plot nerve
plot!(cplx, D, linewidth = 2) # or plot(cplx)
```

### Mapper

Mapper algorithm was described in
"Topological Methods for the Analysis of High Dimensional Data Sets and 3D Object Recognition"
by Gurjeet Singh, Facundo Mémoli and Gunnar Carlsson. Here is an example based
on the description from the Figure 1 of the above paper.

```julia
using TDA, Plots

# generate and plot some dataset
X = hcat(TDA.circlepoints(500, 0.5, noise=0.05)...)'
plot(X[1,:], X[2,:], seriestype=:scatter)

# define Mapper filter function for dataset: f(x) = ||x.x - p.x||
fltfn = (data)->vec(mapslices(p->p[1]-minimum(data[1,:]), data, dims=1))
# plot data colored by filter function values
plot(X[1,:], X[2,:], label="", zcolor=fltfn(X), seriestype=:scatter, ms=2)

# call Mapper algorithm with the particular filter function.
mpr = TDA.mapper(X, filter=fltfn, seed=0, intervals=5, overlap=0.2)

# plot topological layout - mapper graph (by default circular layout is used)
plot(mpr, c=:viridis)
# use `constant_layout` for positioning Mapper graph vertices
# at centers of cover patches
plot(mpr, c=:viridis, complex_layout=TDA.constant_layout)
```

## TODO

- [ ] Plots
    - [x] Persistance Diagram
    - [x] Barcode
    - [x] 1D Simplicial Subcomplex (Graph)
    - [ ] Landscape
- [ ] Mapper
    - [ ] Clustering
        - [x] K-means
        - [ ] Hierarchical
    - [ ] Mode filter functions
    - [x] Plots