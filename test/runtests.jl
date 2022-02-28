################################################################################
#  Copyright 2022, Tom Van Acker (BASF Antwerp)                                #
################################################################################
# TFT.jl                                                                       #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TFT.jl                                    #
################################################################################

using TFT
using Test

f   = 50.0 
ω   = 2 * pi * f
t   = 0.0:0.0001:1.0

@testset "TFT.jl" begin
    
    include("prob/dtft.jl")

end
