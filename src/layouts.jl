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
function spring_layout(cplx::SimplicialComplex; seed::Integer=0, xs=nothing, ys=nothing, kwargs...)
    N = size(cplx, 0)
    rng = Random.MersenneTwister(seed)
    spring_layout(cplx, xs === nothing ? 2 .* rand(rng, N) .- 1.0 : xs,
                        ys === nothing ? 2 .* rand(rng, N) .- 1.0 : ys ; kwargs...)
end

function spring_layout(mpr::Mapper; seed::Integer=0, random::Bool=false, kwargs...)
    xs = random ? nothing : mpr.centers[1,:]
    ys = random ? nothing : mpr.centers[2,:]
    spring_layout(mpr.complex; xs=xs, ys=ys, kwargs...)
end

"""
    circular_layout(cplx)

Arranges verticies of a 1D simplicial complex on a circle.
"""
circular_layout(cplx::SimplicialComplex) = circlepoints(size(cplx, 0)+1, 1.0)
circular_layout(mpr::Mapper) = circular_layout(mpr.complex)

function constant_layout(data::AbstractMatrix)
    inner_layout(cplx::SimplicialComplex) = (data[1,:], data[2,:])
    return inner_layout
end

constant_layout(mpr::Mapper) = (mpr.centers[1,:], mpr.centers[2,:])
