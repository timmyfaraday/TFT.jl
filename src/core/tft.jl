################################################################################
#  Copyright 2022, Tom Van Acker (BASF Antwerp)                                #
################################################################################
# TFT.jl                                                                       #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TFT.jl                                    #
################################################################################

"""
    TFT.tft(prob::AbstractDTFTProblem)

Input:
- `prob::AbstractDTFTProblem`   | DTFT problem struct

Output:
- `sol::AbstractDTFTSolution`   | DTFT solution struct
"""
tft(prob::AbstractDTFTProblem) =
    DTFTSolution(Dict(nh => harmonic_estimator(prob, nh) for nh in prob.h), prob)

"""
    TFT.tft(s::Vector{<:Real}, t::Vector{<:Real}, h::Vector{<:Int}, D::Int, F::Real, K::Int)

Input:
- `s::Vector{<:Real}`   | discrete signal [?]
- `t::Vector{<:Real}`   | discrete time [s]
- `h::Vector{<:Int}`    | harmonic numbers [-]
- `D::Int`              | maximum degree of the derivatives [-]
- `F::Real`             | fundamental frequency [Hz]
- `K::Int`              | degree of the o-spline [-]

Output:
- `sol::DTFTSolution`   | DTFT solution struct
"""
tft(s::Vector{<:Real}, t::Vector{<:Real}, h::Vector{<:Int}, D::Int, F::Real, K::Int) =
    tft(build_problem(s, t, h, D, F, K))