import Clustering

"""
    Mapper result type
"""
struct Mapper
    complex::SimplicialComplex
    filter::Vector{<:Real}
    patches::Vector{Vector{Int}}
    centers::Matrix{<:Real}
end

Base.show(io::IO, mpr::Mapper) = show(io, "Mapper[$(length(mpr.patches))]")

"""
    mapper(X::AbstractMatrix{<:Real}[; filter = eccentricity])

Calculate a simplicial complex of `X` dataset using Mapper algorithm.
"""
function mapper(X::AbstractMatrix{<:Real}; kwargs...)

    # setup parameters
    filterfn = eccentricity
    coverfn = balanced
    clusterfn = Clustering.kmeans
    clusterselectionfn = silhouette
    for (p,v) in kwargs
        if p == :filter
            @assert isa(v, Function) "`$p` parameter must be function"
            filterfn = v
        elseif p == :cover
            @assert isa(v, Function) "`$p` parameter must be function"
            coverfn = v
        elseif p == :clustering
            @assert isa(v, Function) "`$p` parameter must be function"
            clusterfn = v
        elseif p == :clustselection
            @assert isa(v, Function) "`$p` parameter must be function"
            clusterselectionfn = v
        elseif p == :seed
            Random.seed!(v)
        end
    end

    # calculate filter range
    flt = filterfn(X)

    # construct cover of the filter range
    covering = coverfn(flt; kwargs...)

    # using clustering algorithm create patches from cover elements
    patches = Vector{Int}[]
    for c in covering
        println(c)
        if length(c) == 1
            push!(patches, c)
        end
        lbls = clusterselectionfn(clusterfn, view(X, :, c); kwargs...)
        plen = length(unique(lbls))
        println(plen)
        for i in 1:plen
            push!(patches, c[findall(isequal(i), lbls)])
        end
    end

    # combine all patches & determine which has overlaps
    P = length(patches)
    cplx = SimplicialComplex([Simplex(i) for i in 1:P]...)
    for i in 1:P
        println("$i => ", patches[i])
        for j in i+1:P
            overlap = intersect(patches[i], patches[j])
            if length(overlap) > 0
                addsimplex!(cplx, Simplex(i, j))
            end
        end
    end

    # calculate centers of cover patches
    cntrs = hcat((mean(view(X, : ,p), dims=2) for p in patches)...)

    return Mapper(cplx, flt, patches, cntrs)
end

@recipe function f(mpr::Mapper; complex_layout = circular_layout,
                   minvsize = 15, maxvsize = 35)

    xpos, ypos = complex_layout(mpr)

    # set image limits
    xlims --> extrema(xpos) .* 1.2
    ylims --> extrema(ypos) .* 1.2

    # show nerve
    for (i,c) in enumerate(cells(mpr.complex, 1))
        idxs = values(c)
        @series begin
            seriestype := :path
            linewidth --> 2
            linecolor --> :black
            label --> ""
            xpos[idxs], ypos[idxs]
        end
    end

    # calculate vertex attribues
    n = length(mpr.patches)
    xcrd = zeros(n)
    ycrd = zeros(n)
    zcol = zeros(n)
    msize = fill(1,n)
    for (i, p) in enumerate(mpr.patches)
        zcol[i] = mean(mpr.filter[p])
        msize[i] = length(p)
        xcrd[i] = xpos[i]
        ycrd[i] = ypos[i]
    end
    manno = map(string, msize)
    smin, smax = extrema(msize)
    srng = smax-smin
    msize = (maxvsize-minvsize).*(msize .- smin)./srng .+ minvsize

    # show nodes
    @series begin
        seriestype := :scatter
        markersize := msize
        label --> ""
        zcolor := zcol
        series_annotations := manno
        xcrd, ycrd
    end
end
