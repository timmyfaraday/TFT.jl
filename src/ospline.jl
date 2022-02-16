"""
    TFT.sample_ospline(K::Int, N::Int)

Function to obtain samples of a Kth-degree o-spline and its first and second 
derivative.

input:
- K::Int | degree of the o-spline [-]
- N::Int | number of samples of the fundamental cycle [-]
output
- φ⁰::Vector{Real} | samples of the Kth-degree o-spline 
- φ¹::Vector{Real} | samples of the first derivative of the Kth-degree o-spline
- φ²::Vector{Real} | samples of the second derivative of the Kth-degree o-spline
"""
function sample_ospline(K::Int, N::Int)

    # define the normalized time range `rU`
    ΔU  = 1 / N
    rU  = 0.0:ΔU:(1.0-ΔU)
    # define the knots range `rK`
    rK  = [-K:-1..., 1:K...]

    # initialize the Φ-matrices
    Φ⁰, Φ¹, Φ² = ones(N, K+1), ones(N, K+1), ones(N, K+1)

    # fill the Φ-matrices
    for ctr in 0:K
        # one-based counter `cnt` wrt zero-based counter `ctr`
        ctn = ctr + 1
        # base o-spline
        u   = -(K + 1) / 2 + ctr .+ rU
        cf  = [1, -rK[ctn]] ./ -rK[ctn]
        for k in 1:K
            Φ⁰[:, ctn] .*= (u .- rK[ctr + k]) ./ -rK[ctr + k]
            if k > 1
                cf  = _DSP.conv(cf, [1, -rK[ctr + k]] ./ -rK[ctr + k])
            end
        end
        # first derivative by Horner scheme 
        # NB: the Julia Polynomials pkg, takes its coefficients in reverse, 
        # hence the need for the double reserve when determining ψ¹ and ψ².
        ψ¹  = reverse(_POL.derivative(_POL.Polynomial(reverse(cf))).coeffs)
        Φ¹[:, ctn] .*= ψ¹[1]
        for k in 2:K
            Φ¹[:, ctn] .*= u 
            Φ¹[:, ctn] .+= ψ¹[k]
        end
        # second derivative by Horner scheme
        ψ²  = reverse(_POL.derivative(_POL.Polynomial(reverse(ψ¹))).coeffs)
        Φ²[:, ctn] .*= ψ²[1]
        for k in 2:K-1
            Φ²[:, ctn] .*= u 
            Φ²[:, ctn] .+= ψ²[k]
        end
    end

    # return the samples of the o-spline and its derivatives: φ⁰, φ¹, φ²
    return  reshape(Φ⁰, ((K+1) * N, 1)), 
            reshape(Φ¹, ((K+1) * N, 1)), 
            reshape(Φ², ((K+1) * N, 1))

end