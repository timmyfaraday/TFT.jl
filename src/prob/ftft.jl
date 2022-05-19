################################################################################
#  Copyright 2022, Jose Antonio de la O Serna (UANL), Tom Van Acker (BASF)     #
################################################################################
# TFT.jl                                                                       #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TFT.jl                                    #
################################################################################

"""
TFT.ftft(prob::AbstractDTFTProblem)

Input:
- `prob::AbstractDTFTProblem`   | DTFT problem struct

Output:
- `sol::AbstractDTFTSolution`   | DTFT solution struct
"""
ftft(prob::AbstractDTFTProblem) = DTFTSolution(fast_estimator(prob), prob)

"""
TFT.ftft(s::Vector{<:Number}, t::Vector{<:Number}, D::Int, F::Real, K::Int)

Input:
- `s::Vector{<:Number}` | discrete signal [?]
- `t::Vector{<:Number}` | discrete time [s]
- `D::Int`              | maximum degree of the derivatives [-]
- `F::Real`             | fundamental frequency [Hz]
- `K::Int`              | degree of the o-spline [-]

Output:
- `sol::DTFTSolution`   | DTFT solution struct
"""
ftft(s::Vector{<:Number}, t::Vector{<:Number}, D::Int, F::Number, K::Int) =
    ftft(build_problem(s, t, D, F, K)) 

"""
TFT.ftft(s::Vector{<:Number}, t::Vector{<:Number}, h::Vector{<:Int}, D::Int, F::Real, K::Int)
    
Input:
- `s::Vector{<:Number}` | discrete signal [?]
- `t::Vector{<:Number}` | discrete time [s]
- `h::Vector{<:Int}`    | harmonic numbers [-]
- `D::Int`              | maximum degree of the derivatives [-]
- `F::Real`             | fundamental frequency [Hz]
- `K::Int`              | degree of the o-spline [-]
    
Output:
- `sol::DTFTSolution`   | DTFT solution struct
"""
ftft(s::Vector{<:Number}, t::Vector{<:Number}, h::Vector{<:Int}, D::Int, F::Number, K::Int) =
    ftft(build_problem(s, t, h, D, F, K)) 