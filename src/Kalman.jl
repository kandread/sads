module Kalman

function sqrtenkf(A::Matrix, d::Vector, HA::Matrix, E::Matrix)
    m, N = size(A)
    n = size(HA, 1)
    d = reshape(d, n, 1)
    Ap = A .- mean(A, 2)
    S = HA .- mean(HA, 2)
    U0, S0, V0 = svd(S)
    S0p = 1 ./ S0
    S0p[N:end] = 0
    X0 = diagm(S0p) * U0' * E
    U1, S1, V1 = svd(X0)
    X1 = U0 * diagm(S0p)' * U1
    y0 = X1' * (d - mean(HA, 2))
    y2 = pinv(eye(N) + diagm(S1.^2)) * y0
    y3 = X1 * y2
    y4 = S' * y3
    Xma = mean(A, 2) .+ Ap * y4
    X2 = pinv(sqrtm(eye(N) .+ diagm(S1.^2))) * X1' * S
    U2, S2, V2 = svd(X2)
    Xpa = Ap * V2 * sqrtm(eye(N) .- diagm(S2)' * diagm(S2))
    return Xma .+ Xpa
end

export sqrtenkf

end
