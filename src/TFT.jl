################################################################################
#  Copyright 2022, Tom Van Acker (BASF Antwerp)                                #
################################################################################
# TFT.jl                                                                       #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TFT.jl                                    #
################################################################################

module TFT

    # import pkg
    import DSP
    import Polynomials
    import FFTW
    
    # using pkg 
    using Unitful

    # pkg constants 
    const _DSP  = DSP
    const _POL  = Polynomials
    const _UF   = Unitful
    const _FFTW = FFTW 

    # include
    include("prob/dtft.jl")

    include("core/estimator.jl")
    include("core/ospline.jl")
    include("core/tft.jl")

    include("util/util.jl")

    # export
    export  DTFTProblem, DTFTSolution

    export  tft

    export  amplitude, phase, ar_phase, frequency, rocof, phasor, ar_phasor, 
            signal, error
    export  a, ϕ, φ, f, r, ξ, ψ

end
