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

## TODO

- [ ] Graphs
    - [x] Persistance diagram
    - [x] Barcode
    - [x] Simplicial Complex Nerve
    - [ ] Landscape
- [ ] Ripser
- [ ] Mapper