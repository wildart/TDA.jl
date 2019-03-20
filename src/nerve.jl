# Plot a simplicial complex nerve (graph)
@recipe function f(cplx::T, args...) where {T<:SimplicialComplex}
    legend --> :none

    cls = cells(cplx, 1)
    data = length(args) > 0 ? args[1] :
         mapslices(p->p./sqrt(sum(p.^2)), randn(length(cls),2), dims=2)

    for (i,c) in enumerate(cls)
        pts = data[values(c),:]
        @series begin
            seriestype := :path
            linecolor --> :black
            label --> "S$i"
            pts[:,1], pts[:,2]
        end
    end
end

