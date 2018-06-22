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
    y2 = pinv(eye(min(n, N)) + diagm(S1.^2)) * y0
    y3 = X1 * y2
    y4 = S' * y3
    Xma = mean(A, 2) .+ Ap * y4
    # X2 = pinv(sqrtm(eye(N) .+ diagm(S1.^2))) * X1' * S
    # U2, S2, V2 = svd(X2)
    # Xpa = Ap * V2 * sqrtm(eye(N) .- diagm(S2)' * diagm(S2))
    return Xma
end

function letkf(A::Matrix, d::Vector, HA::Matrix, E::Matrix, mapobs)
    rho = 1.05
    lg = length(d)
    mg, k = size(A)
    Yb_avg = mean(HA, 2)
    Xb_avg = mean(A, 2)
    Yb = HA .- Yb_avg
    Xb = A .- Xb_avg
    Aa = zeros(size(A))
    R = (1 / (k - 1)) * (E .- mean(E, 2)) * (E .- mean(E, 2))'
    for i in 1:mg
        j = mapobs(i)
        Rl = zeros(length(j), length(j))
        Rl[:] = R[j, j]
        Ybl = reshape(Yb[j, :], length(j), k)
        C = Ybl' * pinv(Rl)
        Pa = pinv((k-1) / rho * eye(k) + C * Ybl)
        Wa = real.(sqrtm((k-1) * Pa))
        wa = Pa * C * (d[j] .- Yb_avg[j])
        Aa[i, :] = Xb[i, :] .* wa .+ Xb_avg[i]
    end
    return median(Aa, 2)
end

export sqrtenkf, letkf

end
