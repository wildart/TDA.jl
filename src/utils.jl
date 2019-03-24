import Statistics: cov, mean
import LinearAlgebra: inv, det, norm
import MultivariateStats

"""
    circlepoints(r, n)

Generate `n` evenly spaced points on circle of radius `r`
"""
function circlepoints(n, r; noise = 0.0)
    t = LinRange(0.0,2π,n)
    x = r .* cos.(t) .+ rand(n).*noise
    y = r .* sin.(t) .+ rand(n).*noise
    return x, y
end

function gaussian(X::AbstractMatrix{<:Real}; dims=2)
    k = size(X, 1)
    Σ = cov(X, dims=dims)
    Σ⁻¹ = inv(Σ)
    μ = mean(X, dims=dims)
    Xc = X.-μ
    exp.(-0.5*Xc'*Σ⁻¹*Xc)./sqrt(det(Σ)*(2π)^k)
end

function density(X::AbstractMatrix{<:Real}; ε=1.0)
    D = MultivariateStats.pairwise((x,y)->exp(-(1/ε)*norm(x-y)^2.0), X)
    SD = sum(D, dims=2)
    vec(SD./sum(SD))
end

function eccentricity(X::AbstractMatrix{<:Real}; p=2.0)
    N = size(X, 2)
    D = MultivariateStats.pairwise((x,y)->norm(x-y,p), X)
    vec(sum(D, dims=2)./N).^(1/p)
end

"""
    balanced(x[; intervals=10, overlap=0.2])

Balanced cover with `intervals` and `overlap` fraction.
The interval boundaries are distributed so that each patch covers the same fraction of the data set.
"""
function balanced(x::AbstractVector{<:Real}; intervals=10, overlap=0.5)
    xmin, xmax = extrema(x)
    xrng = xmax - xmin
    ilen = 1.0 / (intervals - (intervals-1)*overlap)
    istep = ilen*(1-overlap)

    # compute ranges of intervals with an overlap
    irngs = [let rmin = (i-1)*istep;
                (rmin, (rmin+ilen)) .* xrng .+ xmin
            end
            for i in 1:intervals]

    # generate indexes for cover elements
    patches = map(i->findall(e->i[1]<=e<=i[2], x), irngs)
    return patches #filter!(p->length(p)>0, patches)
end

"""
    adjacency_matrix(cplx)

Construct an adjacency matrix of a one-dimensional simplicial subcomplex of the complex `cplx`.
"""
function adjacency_matrix(cplx::SimplicialComplex)
    N = size(cplx, 0)
    adj = zeros(eltype(celltype(cplx)()),N,N)
    for c in cells(cplx, 1)
        idxs = values(c)
        adj[idxs, idxs] .= 1
    end
    return adj
end

"""

    spring_layout(cplx[; C=2.0, maxiter=100, seed=0])

Calculate a 1D simplicial complex layout using spring/repulsion model.

This borrows with modifications from [GraphPlot](https://github.com/JuliaGraphs/GraphPlot.jl)

# Arguments
- `cplx::SimplicialComplex`: the simplicial complex.
- `C::Real`: the constant to fiddle with density of resulting layout.
- `maxiter::Integer`: the number of iterations we apply the forces.
- `starttemp::Real`: the initial "temperature", controls movement per iteration.
- `seed::Integer`: the seed for pseudorandom generation of locations (default = 0).
"""
function spring_layout(cplx::SimplicialComplex, locs_x, locs_y; C=2.0, maxiter=100, starttemp=2.0) where {T<:Integer}

    N = size(cplx, 0)
    adj_matrix = adjacency_matrix(cplx)

    # The optimal distance bewteen vertices
    K = C * sqrt(4.0 / N)

    # Store forces and apply at end of iteration all at once
    force_x = zeros(N)
    force_y = zeros(N)

    @inbounds for iter = 1:maxiter
        # Calculate forces
        for i = 1:N
            force_vec_x = 0.0
            force_vec_y = 0.0
            for j = 1:N
                i == j && continue
                d_x = locs_x[j] - locs_x[i]
                d_y = locs_y[j] - locs_y[i]
                d   = sqrt(d_x^2 + d_y^2)
                if adj_matrix[i,j] != zero(eltype(adj_matrix)) || adj_matrix[j,i] != zero(eltype(adj_matrix))
                    # F = d^2 / K - K^2 / d
                    F_d = d / K - K^2 / d^2
                else
                    # Just repulsive
                    # F = -K^2 / d^
                    F_d = -K^2 / d^2
                end
                # d  /          sin θ = d_y/d = fy/F
                # F /| dy fy    -> fy = F*d_y/d
                #  / |          cos θ = d_x/d = fx/F
                # /---          -> fx = F*d_x/d
                # dx fx
                force_vec_x += F_d*d_x
                force_vec_y += F_d*d_y
            end
            force_x[i] = force_vec_x
            force_y[i] = force_vec_y
        end
        # Cool down
        TEMP = starttemp / iter
        # Now apply them, but limit to temperature
        for i = 1:N
            force_mag  = sqrt(force_x[i]^2 + force_y[i]^2)
            scale      = min(force_mag, TEMP)/force_mag
            locs_x[i] += force_x[i] * scale
            #locs_x[i]  = max(-1.0, min(locs_x[i], +1.0))
            locs_y[i] += force_y[i] * scale
            #locs_y[i]  = max(-1.0, min(locs_y[i], +1.0))
        end
    end

    # Scale to unit square
    min_x, max_x = minimum(locs_x), maximum(locs_x)
    min_y, max_y = minimum(locs_y), maximum(locs_y)
    function scaler(z, a, b)
        2.0*((z - a)/(b - a)) - 1.0
    end
    map!(z -> scaler(z, min_x, max_x), locs_x, locs_x)
    map!(z -> scaler(z, min_y, max_y), locs_y, locs_y)

    return locs_x, locs_y
end

import Random
function spring_layout(cplx::SimplicialComplex; seed::Integer=0, kws...)
    N = size(cplx, 0)
    rng = Random.MersenneTwister(seed)
    spring_layout(cplx, 2 .* rand(rng, N) .- 1.0, 2 .* rand(rng,N) .- 1.0; kws...)
end

"""
    circular_layout(cplx)

Arranges verticies of a 1D simplicial complex on a circle.
"""
circular_layout(cplx::SimplicialComplex) = circlepoints(size(cplx, 0)+1, 1.0)

function constant_layout(data::AbstractMatrix)
    inner_layout(cplx::SimplicialComplex) = (data[1,:], data[2,:])
    return inner_layout
end
