"""
    TFT.sample_ospline(K::Int, N::Int)
`sample_ospline` is a function to obtain the samples the `K`th O-Spline for `N` samples per cycle
"""
function sample_ospline(K::Int, N::Int)

    # define the normalized time range `rU`
    Δu  = 1 / N
    rU  = 0.0:Δu:(1.0-Δu)
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
        Ψ¹  = _POL.derivative(_POL.Polynomial(cf)).coeffs
        Φ¹[:, ctn] .*= Ψ¹[1]
        for k in 2:K
            Φ¹[:, ctn] .*= u 
            Φ¹[:, ctn] .+= Ψ¹[k]
        end
        # second derivative by Horner scheme
        Ψ²  = _POL.derivative(_POL.Polynomial(Ψ¹)).coeffs
        Φ²[:, ctn] .*= Ψ²[1]
        for k in 2:K-1
            Φ²[:, ctn] .*= u 
            Φ²[:, ctn] .+= Ψ²[k]
        end
    end

    # return the ospline sample 
    return reshape(Φ⁰, ((K+1) * N, 1)), reshape(Φ¹, ((K+1) * N, 1)), reshape(Φ², ((K+1) * N, 1))

end