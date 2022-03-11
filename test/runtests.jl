################################################################################
#  Copyright 2022, Tom Van Acker (BASF Antwerp)                                #
################################################################################
# TFT.jl                                                                       #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TFT.jl                                    #
################################################################################

using TFT
using Test

# fundamental frequency and angular frequency
F   = 50.0 
ω   = 2 * pi * F

# tft input
D   = 2
K   = 9

# discrete time
t   = 0.0:0.0001:1.0
tm  = 0.5
idm = findfirst(x -> x == 0.5, t)

# tolerances
atol = 1e-6

@testset "TFT.jl" begin
    
    include("prob/dtft.jl")

end
