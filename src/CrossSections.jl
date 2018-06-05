module CrossSections

abstract type CrossSection end

mutable struct Rectangular <: CrossSection 
    b::Real
    y::Real
    S0::Real
    n::Real
end

function area(xs::Rectangular)
    return xs.b * xs.y
end

function perimeter(xs::Rectangular)
    return 2 * xs.y + xs.b
end

function topwidth(xs::Rectangular)
    return xs.b
end

mutable struct Trapezoidal <: CrossSection
    b::Real
    y::Real
    m::Real
    S0::Real
    n::Real
end

function area(xs::Trapezoidal)
    return xs.y * (xs.b + xs.m * xs.y)
end

function perimeter(xs::Trapezoidal)
    return xs.b + 2 * xs.y * (1 + xs.m^2)^0.5
end

function topwidth(xs::Trapezoidal)
    return xs.b + 2 * xs.m * xs.y
end

mutable struct Triangular <: CrossSection
    y::Real
    m::Real
end

function area(xs::Triangular)
    return xs.m * xs.y^2
end

function perimeter(xs::Triangular)
    return 2 * xs.y * (1 + xs.m^2)^0.5
end

function topwidth(xs::Triangular)
    return 2 * xs.m * xs.y
end

mutable struct Parabolic <: CrossSection
    b1::Real
    y::Real
    y1::Real
end

function area(xs::Parabolic)
    return (2/3) * topwidth(xs) * xs.y
end

function perimeter(xs::Parabolic)
    B = topwidth(xs)
    x = 4 * xs.y / B
    return (B / 2) * (1 + x^2)^0.5 + (1 / x) * log(x + (1 + x^2)^0.5)
end

function topwidth(xs::Parabolic)
    return xs.b1 * (xs.y / xs.y1)^0.5
end

function depth(xs::CrossSection)
    return area(xs) / topwidth(xs)
end

function radius(xs::CrossSection)
    return area(xs) / perimeter(xs)
end

export Rectangular, Trapezoidal, Triangular, Parabolic

export area, radius, depth, topwidth, perimeter

end
