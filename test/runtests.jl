################################################################################
#  Copyright 2022, Tom Van Acker (BASF)                                        #
################################################################################
# TFT.jl                                                                       #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TFT.jl                                    #
################################################################################

# using pkgs
using ForwardDiff
using Test
using TFT

# pkg constants
const _FD = ForwardDiff

# tolerances
atol = 1e-6

@testset "TFT.jl" begin

    include("prob/tft.jl")
    include("prob/ftft.jl")

end
