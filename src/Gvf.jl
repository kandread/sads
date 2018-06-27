module Gvf

using CrossSections
using DifferentialEquations

"""Calculate the Froude number."""
function froude(Q::Real, xs::CrossSection)
    g = 9.806
    Fr = Q / (area(xs) * sqrt(g * depth(xs)))
    return Fr
end

"""Ordinary differential equation describing the Gradually-Varied-Flow model."""
function dydx(y, p, x)
    i = trunc(Int, ceil(x / (p[3] / length(p[1]))))
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
function gvf(Q, ybc, S0, n, x, b)
    c = [Rectangular(b[i], ybc, S0[i], n) for i in 1:length(x)]
    # c = [Trapezoidal(b[i], ybc, 1.0, S0[i], n) for i in 1:length(x)]
    # prob = ODEProblem(dydx, ybc, (x[1], x[end]), (Q, c))
    prob = DiscreteProblem(dydx, ybc, (x[1], x[end]), (Q, c, x[end]))
    h = try
        sol = solve(prob, Tsit5(), abstol=1e-2, saveat=x)
        sol.u
    catch
        zeros(length(x)) - 9999.
    end
    return h
end

"""Solve Gradually-Varied-Flow equation using the second-order
predictor-corrector method."""
function cp2(Q, ybc, S0, n, x, b)
    abstol = 1e-2
    y = zeros(length(x))
    y[1] = ybc
    c = [Rectangular(b[i], ybc, S0[i], n) for i in 1:length(x)]
    f = dydx(ybc, (Q, c, x[end]), x[1])
    for j in 2:length(x)
        init = true
        f1 = f
        while (abs(f1 - f) > abstol || init)
            init = false
            f = f1
            y[j] = y[j-1] + (f + f1) / 2 * (x[j] - x[j-1])
            f1 = dydx(y[j], (Q, c, x[end]), x[j])
        end
    end
    return y
end

export gvf

end
