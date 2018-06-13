module Gvf

using CrossSections

"""Calculate the Froude number."""
function froude(Q::Real, xs::CrossSection)
    g = 9.806
    Fr = Q / (area(xs) * sqrt(g * depth(xs)))
    return Fr
end

"""Ordinary differential equation describing the Gradually-Varied-Flow model."""
function dydx(y, p, x)
    i = trunc(Int, ceil(x / (6000 / length(p[1]))))
    i = i > 0 ? i : 1
    # FIXME: Change this dynamically depending on whether the flow is subcritical or supercritical
    pm = -1  # integrate upstream
    Q = p[1]
    xs = p[2][i]
    xs.y = y
    S0 = xs.S0
    Sf = xs.n^2 * Q^2 / (area(xs)^2 * radius(xs)^(4/3))
    Fr = froude(Q, xs)
    pm * (S0 - Sf) / (1 - Fr^2)
end

"""
Calculate a water surface profile by solving the Gradually-Varied-Flow equation.
`Q`: discharge
`ybc`: downstream boundary condition for depth
`z0`: bed elevation at downstream cross section
`S0`: bed slope for reach
`n`: Manning's roughness coefficient
`x`: downstream distance for each cross section
`b`: bottom width for each cross section
"""
function gvf(Q, ybc, z0, S0, n, x, b)
    c = [Rectangular(b[i], ybc, S0, n) for i in 1:length(x)]
    prob = ODEProblem(dydx, ybc, (x[1], x[end]), (Q, c))
    sol = solve(prob, Tsit5(), abstol=1e-2, saveat=x)
    h = sol.u[end:-1:1]
    return h
end

export gvf

end
