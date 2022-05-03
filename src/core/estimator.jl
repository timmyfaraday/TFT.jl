################################################################################
#  Copyright 2022, Tom Van Acker (BASF Antwerp)                                #
################################################################################
# TFT.jl                                                                       #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TFT.jl                                    #
################################################################################

"""
    TFT.harmonic_estimator(prob::AbstractDTFTProblem, H::Int)

Function to obtain the up-to-Dth-degree derivative of the Hth-harmonic dynamic 
phasor.

Input:
- `prob::AbstractDTFTProblem`   | DTFT problem struct
- `H::Int`                      | harmonic number [-]

Output:
- `X::Matrix{<:Complex}`        | up-to-Dth-degree derivative of the 
                                | Hth-harmonic dynamic phasor
"""
function harmonic_estimator(prob::AbstractDTFTProblem, H::Int)
    # define the extended normalized time range `rU`
    Δu  = 1 / prob.N
    rU  = 0.0:Δu:(prob.K + 1 - Δu)

    # samples of up-to-Dth-degree derivative of the Kth-degree o-spline `Φ`
    Φ   = sample_ospline(prob.D, prob.K, prob.N)
    # update `Φ` to reflect the samples of up-to-Dth-degree derivative of the 
    # Hth-harmonic bandpass filter
    Y   = Φ .* exp.((2 * pi * im * H) .* rU) ./ prob.N .* 
            _UF.ustrip(prob.F).^collect(0:prob.D)'

    # define the appropriate range in the convolution `rC` using an offset `O`:
    # if (K+1)*N is even    → O = ((K+1)*N) / 2 + 1
    # if (K+1)*N is oneven  → O = ((K+1)*N - 1) / 2 + 1
    # NB: the Julia DSP pkg does not allow for the keyword 'same' in its conv-
    # function, as is the case in MATLAB, the range `rC` mimics that behavior.
    O   = floor(Int, ((prob.K + 1) * prob.N) / 2) + 1 
    rC  = O:O+length(prob.s)-1
    
    # return the up-to-Dth-degree derivative of the H-th harmonic dynamic phasor
    uS = _UF.unit(first(prob.s))
    return _DSP.conv(_UF.ustrip.(prob.s), Y)[rC, :]uS
end