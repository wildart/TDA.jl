import Clustering

"""
    Mapper result type
"""
struct Mapper
    complex::SimplicialComplex
    filter::Vector{<:Real}
    patches::Vector{Vector{Int}}
end

Base.show(io::IO, mpr::Mapper) = show(io, "Mapper[$(length(mpr.patches))]")

"""
    centers(mpr, X)

Return centers of cover patches for a Mapper complex `mpr` calculated from dataset `X`.
"""
function centers(mpr::Mapper, X::AbstractMatrix{<:Real})
    return [mean(view(X, : ,p), dims=2) for p in mpr.patches]
end

"""
    mapper(X::AbstractMatrix{<:Real}[; filter = eccentricity])

Calculate a simplicial complex of `X` dataset using Mapper algorithm.
"""
function mapper(X::AbstractMatrix{<:Real};
                filter::Function = eccentricity,
                cover::Function = balanced,
                clustering::Function = Clustering.kmeans,
                cluster_selection::Function = silhouette)

    # calculate filter range
    flt = filter(X)

    # construct cover of the filter range
    covering = cover(flt)

    # using clustering algorithm create patches from cover elements
    patches = Vector{Int}[]
    for c in covering
        lbls = cluster_selection(clustering, view(X, :, c))
        for i in 1:length(unique(lbls))
            push!(patches, c[findall(isequal(i), lbls)])
        end
    end

    # combine all patches & determine which has overlaps
    P = length(patches)
    cplx = SimplicialComplex([Simplex(i) for i in 1:P]...)
    for i in 1:P
        for j in i+1:P
            overlap = intersect(patches[i], patches[j])
            if length(overlap) > 0
                addsimplex!(cplx, Simplex(i, j))
            end
        end
    end

    return Mapper(cplx, flt, patches)
end

@recipe function f(mpr::Mapper; complex_layout=circular_layout)

    xpos, ypos = complex_layout(mpr.complex)

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

    # show nodes
    for (i, p) in enumerate(mpr.patches)
        mcolor = mean(mpr.filter[p])
        msize = length(p)
        @series begin
            seriestype := :scatter
            markersize := msize
            label --> "$msize"
            zcolor := mcolor
            [xpos[i]], [ypos[i]]
        end
    end
end
