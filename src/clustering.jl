import Clustering
import Distances
import LinearAlgebra: dot

"""
    noselection(algo, data)

Perform no clustering selection.
"""
noselection(algo::Function, data::AbstractMatrix{<:Real}) = Clustering.assignments(algo(data))

"""
Knee point detection
"""
function knee(x)
    n = length(x)
    pts = hcat(1:n, x)'
    # calculate line vector
    lvec = pts[:,end] - pts[:,1]
    lvec ./= sqrt(sum(lvec.^2))
    # find a vector parallel to the line vector
    vecfst = pts .- pts[:,1]
    vecfstpar = lvec * mapslices(v->dot(v, lvec), vecfst, dims=1)
    # calculate a length of a projection to the line vector
    ldist = vec(sqrt.(sum((vecfst - vecfstpar).^2, dims=1)))
    # get maximal projection
    return findmax(ldist)
end

function getmaxk(N::Int; kwargs...)
    N <= 5 && return 2

    # setup parameters
    Kmax = 10 # default
    for (p,v) in kwargs
        if p == :maxk
            Kmax = v
        end
    end
    return min(floor(Int, N/2), Kmax)
end

"""
    elbow(algo, data)

Perform automatic selection of the best clustering of `data` by `algo` using elbow method.
Clustering algorithm `algo` must be a function with one argument which determines number of clusters.
"""
function elbow(algo::Function, data::AbstractMatrix{<:Real}; kwargs...)
    N = size(data,2)
    Kmax = getmaxk(N; kwargs...)

    S = [let C = algo(data, k)
            # println("$k => ", C.centers)
            # println("$k => ", C.costs)
            C.totalcost => Clustering.assignments(C)
         end
    for k in 2:Kmax]
    pushfirst!(S, sum(mapslices(p->sum(p.^2), data .- mean(data, dims=2), dims=1))=>fill(1,N))

    # return a clustering corespondig to an elbow point
    return S[knee(map(first, S)) |> last] |> last
end

"""
    silhouette(algo, data)

Perform automatic selection of the best clustering of `data` by `algo` using silhouette method.
Clustering algorithm `algo` must be a function with one argument which determines number of clusters.
"""
function silhouette(algo::Function, data::AbstractMatrix{<:Real}; kwargs...)
    N = size(data,2)
    Kmax = getmaxk(N; kwargs...)
    S = [let C = algo(data, k),
             D = Distances.pairwise(Distances.Euclidean(), data)
             mean(Clustering.silhouettes(C, D)) => Clustering.assignments(C)
    end
    for k in 2:Kmax]

    # return a clustering corespondig to a maximum average silhouette
    return S[findmax(map(first, S)) |> last] |> last
end
