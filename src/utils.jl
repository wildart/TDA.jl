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
function balanced(x::AbstractVector{<:Real}; kwargs...)

    # setup parameters
    intervals=10
    overlap=0.5
    for (p,v) in kwargs
        if p == :intervals
            intervals = v
        elseif p == :overlap
            overlap = v
        end
    end

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
    return patches
end
