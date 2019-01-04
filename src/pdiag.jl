# Plot intervals as a persistance diagram or barcode
@recipe function f(ints::Vector{Vector{T}}; skipzero = false) where {T<:Interval}
    seriestype --> :diagram

    maxval = mapreduce(l->maximum(filter(!isinf, map(i->i.d,l))), max, ints)
    minval = mapreduce(l->minimum(filter(!isinf, map(i->i.b,l))), min, ints)
    #col = vcat(map(l->map(i->isinf(i.d) ? :red : :auto ,l), ints)...)

    if plotattributes[:seriestype] == :barcode # Persistance Barcode
        xticks := minval:maxval
        xlims  := (minval, maxval+0.05)
        yticks := nothing
        ylims := (0.0, 1.05)
        legend := :none
        layout := (length(ints),1)
        for d in 1:length(ints)
            points = map(i -> i.b == i.d, ints[d])
            step = 1 / (skipzero ? sum(.!points) : length(ints[d]))
            y = step
            for (i, ispoint) in zip(ints[d], points)
                xs, ys = [i.b, isinf(i.d) ? maxval : i.d], [y, y]
                println("ispoint: $ispoint")
                println(!ispoint && skipzero)
                ispoint && skipzero && continue
                @series begin
                    title := "Degree $(d-1)"
                    linewidth --> 2
                    subplot := d
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
        xs = map(l->map(i->i.b,l), ints)
        ys = map(l->map(i->isinf(i.d) ? maxval : i.d,l), ints)
        mkr = vcat(map(l->map(i->isinf(i.d) ? :rect : :circle,l), ints)...)
        legend --> :bottomright
        legendtitle --> "Homology"
        label --> ["Degree $(i-1)" for i in 1:length(ints)]
        xticks := minval:maxval
        xlims := (minval-0.3, maxval+0.3)
        xlabel := "birth"
        yticks := 0:maxval
        ylabel := "death"
        ylims := (minval, maxval+0.3)
        @series begin
            primary := false
            seriestype := :path
            linecolor := :black
            [0, maxval], [0, maxval]
        end

        seriestype := :scatter
        markershape := mkr
        #markerstrokecolor := col
        xs, ys
    end
end
