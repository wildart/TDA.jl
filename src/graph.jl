# Plot a 1D simplicial subcomplex (graph) of `cplx`
@recipe function f(cplx::T, args...) where {T<:SimplicialComplex}
    data = length(args) > 0 ? args[1] :
                              circlepoints(size(cplx, 0)+1, 1.0)
    for (i,c) in enumerate(cells(cplx, 1))
        pts = data[values(c),:]
        @series begin
            seriestype := :path
            linewidth --> 2
            linecolor --> :black
            label --> "S$i"
            pts[:,1], pts[:,2]
        end
    end
end
