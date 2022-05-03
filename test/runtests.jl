################################################################################
#  Copyright 2022, Tom Van Acker (BASF Antwerp)                                #
################################################################################
# TFT.jl                                                                       #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TFT.jl                                    #
################################################################################

# using pkgs
using ForwardDiff
using Test
using TFT
using Unitful

# pkg constants
const _FD = ForwardDiff

# tolerances
atol = 1e-6

@testset "TFT.jl" begin
    
    include("unitful.jl")

    include("prob/dtft.jl")

end
