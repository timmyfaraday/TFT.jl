# use pkgs
using DSP, Polynomials

# input
## fundamental frequency [Hz]
F₁ = 50.0
## fundamental period [s]
T₁ = 1 / F₁
## range of harmonics 
H = 0:1:50
## number of samples per fundamental period [-]
N = 100
## degree of the Taylor polynomial [-]
K = 2

function ospline(K::Int, N::Int)
    Δ       = 1 / N
    unit    = 0:Δ:1-Δ
    knots   = -K:K

    P₀, P₁, P₂   = ones(N, K+1), ones(N, K+1), ones(N, K+1)

    ir0 = 0 
    for hnint in 1:K+1
        u = -(K+1) / 2 - 1 hnint + unit
        coef = [1, -knots[ir0+1]] / -knots[ir0+1]
        for k in 1:K
            P₀[:,hnint] = P₀[:,hnint] .* (u - knots[ir0+k]) / (-knots[ir0+k])
            if k > 1
                coef = DSP.conv(coef, [1, -knots[ir0+k]] / (-knots[ir0+k]))
            end
        end

        C₀, C₁, C₂ = zeros(K+1, hnint), zeros(K, hnint), zeros(K-1, hnint) 
        
        C₀[:,hnint] = coef
        C₁[:,hnint] = polyder(Poly(coef))
        C₂[:,hnint] = polyder(Poly(C₁[:,hnint]))

        # first derivative by Horner scheme
        P₁[:,hnint] = C₁[1,hnint] * P₁[:,hnint]
        for k in 1:K-1
            P₁[:,hnint] = P₁[:,hnint] .* u + C₁[k+1,hnint]
        end

        # second derivative by Horner scheme
        P₂[:,hnint] = C₂[1,hnint] * P₂[:,hnint]
        for k in 1:K-2
            P₂[:,hnint] = P₂[:,hnint] .* u + C₂[k+1,hnint]
        end

        ir0 = ir0 + 1
    end

    Φ₁ = reshape(P₁, [(K+1)*N, 1])
    Φ₂ = reshape(P₂, [(K+1)*N, 1])
    Φ₃ = reshape(P₃, [(K+1)*N, 1])

    return Φ₁, Φ₂, Φ₃
end