"""
    circlepoints(r, n)

Generate `n` evenly spaced points on circle of radius `r`
"""
function circlepoints(n, r)
    t = LinRange(0.0,2π,n)
    x = r .* cos.(t)
    y = r .* sin.(t)
    hcat(x, y)
end
