module Sads

using DifferentialEquations
using CrossSections
using Kalman

function froude(Q::Real, xs::CrossSection)
    g = 9.806
    Fr = Q / (area(xs) * sqrt(g * depth(xs)))
    return Fr
end

function dydx(y, p, x)
    i = trunc(Int, ceil(x / (6000 / length(p[1]))))
    i = i > 0 ? i : 1
    pm = -1  # integrate upstream
    Q = p[1]
    xs = p[2][i]
    xs.y = y
    S0 = xs.S0
    Sf = xs.n^2 * Q^2 / (area(xs)^2 * radius(xs)^(4/3))
    Fr = froude(Q, xs)
    pm * (S0 - Sf) / (1 - Fr^2)
end

function gvf()
    
end

function estimate()
    
end

end
