# Plot intervals as a persistance diagram or barcode
@recipe function f(ints::Dict{Int, Vector{Interval}};
        maxoutdim=length(ints),
        skipzero = false)
    # set default plot type
    seriestype --> :diagram

    vals = vcat(([i.b, i.d] for i in vcat(values(ints)...))...)
    minval, maxval = extrema(filter!(!isinf, vals))
    #col = vcat(map(l->map(i->isinf(i.d) ? :red : :auto ,l), ints)...)
    dims = [d for d in sort!(filter!(d->d<=maxoutdim, collect(keys(ints))))]

    if plotattributes[:seriestype] == :barcode # Persistance Barcode
        #xticks := minval:maxval
        xlims  := (minval, maxval*1.01)
        yticks := nothing
        ylims := (0.0, 1.05)
        legend := :none
        layout := (length(dims),1)
        for (pi, d) in enumerate(dims)
            dmaxval = max(maxval, foldl(max, filter(!isinf, map(i->i.d, ints[d])); init=-Inf))
            dminval = min(minval, foldl(min, filter(!isinf, map(i->i.b, ints[d])); init=Inf))
            points = map(i -> i.b == i.d, ints[d])
            step = 1 / (skipzero ? sum(.!points) : length(ints[d]))
            y = step
            for (i, ispoint) in zip(ints[d], points)
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
        xs = map(l->map(i->i.b,l), ints[d] for d in dims)
        ys = map(l->map(i->isinf(i.d) ? maxval : i.d,l), ints[d] for d in dims)
        mkr = map(l->map(i->isinf(i.d) ? :rect : :circle,l), ints[d] for d in dims)
        legend --> :bottomright
        legendtitle --> "Homology"
        label --> ["Degree $(i-1)" for i in 1:length(ints)]
        #xticks := minval:maxval
        xlims := (minval-padding, maxval+padding)
        xlabel := "birth"
        #yticks := 0:maxval
        ylabel := "death"
        ylims := (minval, maxval*1.03)

        @series begin
            primary := false
            seriestype := :path
            linecolor := :black
            [0, maxval], [0, maxval]
        end

        for i in 1:length(xs)
            idxs = skipzero ? xs[i] .!== ys[i] : .!BitArray(undef, length(xs[i]))
            @series begin
                primary := true
                seriestype := :scatter
                markershape := mkr[i][idxs]
                markersize --> 2
                #markerstrokecolor := col
                xs[i][idxs], ys[i][idxs]
            end
        end
    end
end
