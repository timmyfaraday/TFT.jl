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
    
    # return the up-to-Dth-degree derivative of the H-th harmonic complex envelope
    return _DSP.conv(prob.s, Y)[rC, :]  
end
"""
TFT.FTFT(Φ::Matrix{<:Real}, D::Int, K::Int, N::Int,s::Matrix{<:Real})

Function to obtain the Fast Taylor Fourier transform FTFT
Uses the samples of the up-to-Dth-degree derivatives of Kth-degree o-spline in matrix Φ.

Input:
- `Φ::Matrix{<:Real}` | Contains inthe rows the O-spline of Kth order and its derivatives up to order D
- `D::Int`  | maximum degree of the derivative [-]
- `K::Int`  | degree of the o-spline [-]
- `N::Int`  | number of samples of the fundamental cycle [-]
- `s::Matrix{<:Real}` | (K+1)Nx1 veector with signal samples  

Output
- `ξ::Matrix{<:Complex}` | (D+1)Nx1 vector with the harmonic Taylor-Fourier coefficients up to the D derivative.
                         | obtained with the Kth-degree o-spline 
"""
function FTFT(Φ::Matrix{<:Real}, D::Int, K::Int, N::Int, s::Matrix{<:Real})
    for nd=1:D+1
        hd=Φ[nd,:]'.*s
        S=reshape(hd,(N,K+1))
        shd=sum(S,dims=2)
        ξ[(nd-1)*N+1:nd*N]=fft(shd)/N
    end    

