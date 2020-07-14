@recipe function f(ints::Vector{T}) where {T<:AbstractInterval}
    Dict(0=>ints)
end

# Plot intervals as a persistance diagram or barcode
@recipe function f(dgm::Dict{Int, Vector{T}};
                   maxoutdim=1,
                   skipzero = false) where {T<:AbstractInterval}
    # set default plot type
    seriestype --> :diagram

    # get interval dimensions
    vals = vcat(([i.b, i.d] for i in vcat(values(dgm)...))...)
    dims = [d for d in sort!(filter!(d->d<=maxoutdim, collect(keys(dgm))))]
    # compute range of intervals
    minval, maxval = extrema(filter!(!isinf, vals))

    if plotattributes[:seriestype] == :barcode # Barcode
        #xticks := minval:maxval
        xlims  := (minval, maxval*1.01)
        yticks := nothing
        ylims := (0.0, 1.05)
        legend := :none
        layout := (length(dims),1)
        for (pi, d) in enumerate(dims)
            dmaxval = max(maxval, foldl(max, filter(!isinf, map(i->i.d, dgm[d])); init=-Inf))
            dminval = min(minval, foldl(min, filter(!isinf, map(i->i.b, dgm[d])); init=Inf))
            points = map(i -> i.b == i.d, dgm[d])
            step = 1 / (skipzero ? sum(.!points) : length(dgm[d]))
            y = step
            idxs = sortperm(dgm[d])
            for (i, ispoint) in zip(dgm[d][idxs], points)
                xs, ys = [i.b, isinf(i.d) ? maxval : i.d], [y, y]
                ispoint && skipzero && continue
                @series begin
                    title := "Degree $d"
                    titlefontsize --> 12
                    linewidth --> 2
                    subplot := pi
                    primary := true
                    arrow := isinf(i.d)
                    seriescolor --> :blue
                    seriestype := ispoint ? :scatter : :path
                    markersize --> 2
                    xs, ys
                end
                y += step
            end
        end
    else # Persistance Diagram
        padding = (maxval - minval)*0.01 # 3%
        xs = map(l->map(i->i.b,l), dgm[d] for d in dims)
        ys = map(l->map(i->isinf(i.d) ? maxval : i.d,l), dgm[d] for d in dims)
        mkr = map(l->map(i->isinf(i.d) ? :rect : :circle,l), dgm[d] for d in dims)
        legend --> :bottomright
        legendtitle --> "Homology"
        # label --> ["Degree $(i-1)" for i in 1:length(dgm)]
        # xticks := minval:maxval
        xlims := (minval-padding, maxval+padding)
        xguide := "birth"
        # yticks := 0:maxval
        yguide := "death"
        ylims := (minval, maxval*1.03)

        @series begin
            primary := false
            seriestype := :path
            linecolor := :black
            [0, maxval], [0, maxval]
        end

        for i in 1:length(xs)
            idxs = skipzero ? xs[i] .!== ys[i] : BitArray(true for i in 1:length(xs[i]))
            @series begin
                primary := true
                seriestype := :scatter
                markershape := mkr[i][idxs]
                markersize --> 2
                label --> "Degree $(i-1)"
                # markerstrokecolor := col
                xs[i][idxs], ys[i][idxs]
            end
        end
    end
end
