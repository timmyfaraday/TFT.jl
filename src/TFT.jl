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

    # pkg constants 
    const _DSP = DSP
    const _POL = Polynomials

    # include
    include("prob/dtft.jl")

    include("core/estimator.jl")
    include("core/ospline.jl")
    include("core/tft.jl")

    include("util/util.jl")

    # export
    export DTFTProblem, DTFTSolution

    export tft

    export amplitude, angle, frequency, rocof, phasor, signal, overall_signal, error

end
