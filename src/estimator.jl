"""
    TFT.harmonic_state_estimator(s::Vector{<:Real}, K::Int, N::Int, h::Int, F::Real)

Function to obtain the hth-harmonic amplitude and angle of a signal `s`, its 
frequency and ROCOF (Rate Of Change Of Frequency).

input:
- s::Vector{<:Real} | samples of the signal
- K::Int | degree of the o-spline [-]
- N::Int | number of samples of the fundamental cycle [-]
- H::Int | harmonic number [-]
- F::Real | fundamental frequency of the signal [Hz]
output:
- a⁰::Vector{Real} | amplitude of the signal
- ϕ⁰::Vector{Real} | angle of the signal
- a¹::Vector{Real} | amplitude of the frequency 
- ϕ¹::Vector{Real} | angle of the frequency 
- a²::Vector{Real} | amplitude of the ROCOF 
- ϕ²::Vector{Real} | angle of the ROCOF
"""
function harmonic_state_estimator(s::Vector{<:Real}, K::Int, N::Int, H::Int, F::Real)

    # define the zero-based index range of the sample `rS`
    S   = length(s)
    rS  = 0:S-1

    # define the extended normalized time range `rU`
    Δu  = 1 / N
    rU  = 0.0:Δu:(K + 1 - Δu)

    # samples of the Kth-degree o-spline and its derivative: `φ⁰`, `φ¹`, `φ²`
    φ⁰, φ¹, φ² = sample_ospline(K, N)
    # samples of the hth-harmonic bandpass filters and its derivatives: `h`, `f`, `r`
    h   = φ⁰ .* exp.((2 * pi * im * H) .* rU) ./ N 
    f   = φ¹ .* exp.((2 * pi * im * H) .* rU) ./ N .* F
    r   = φ² .* exp.((2 * pi * im * H) .* rU) ./ N .* F^2

    # define the appropriate range in the convolution `rC` using an offset `O`:
    # if (K+1)*N is even    → O = ((K+1)*N) / 2 + 1
    # if (K+1)*N is oneven  → O = ((K+1)*N - 1) / 2 + 1
    # NB: the Julia DSP pkg does not allow for the keyword same in its conv-
    # function, as is the case in MATLAB, the range `rC` mimics that 
    # functionality.
    O   = floor(Int, ((K + 1) * N)) / 2 + 1 
    rC  = O:O+S-1
    
    # amplitude and phase estimation of the signal
    ξ⁰  = _DSP.conv(s, h)[rC]
    ah⁰ = abs.(ξ⁰)
    ϕh⁰ = angle.(ξ⁰)
    # amplitude and phase estimation of the frequency
    ξ¹  = _DSP.conv(s, f)[rC] .* exp.(-im .* ϕh⁰)
    ah¹ = real.(ξ¹)
    ϕh¹ = imag.(ξ¹) ./ ah⁰
    # amplitude and phase estimation of the ROCOF
    ξ²  = _DSP.conv(s, r)[rC] .* exp.(-im * ϕh⁰)
    ah² = real.(ξ²) .+ ah⁰ .* ϕh¹.^2
    ϕh² = (imag.(ξ²) .- 2 .* ah¹ .* ϕh¹) ./ ah⁰    

    # return a⁰, ϕ⁰, a¹, ϕ¹, a², ϕ²
    return  2 .* ah⁰, angle.(ξ⁰ .* exp.((-2 * pi * im / N) .* rS)), 
            ah¹, ϕh¹ ./ (2 * pi), 
            ah², ϕh² ./ (2 * pi)

end