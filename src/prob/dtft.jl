################################################################################
#  Copyright 2022, Tom Van Acker (BASF Antwerp)                                #
################################################################################
# TFT.jl                                                                       #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TFT.jl                                    #
################################################################################

# abstract problem types
"""
    TFT.AbstractDTFTProblem
    
Base for types which define discrete taylor-fourier transform problems.
"""
abstract type AbstractDTFTProblem end

# problem types
"""
    TFT.DTFTProblem

Defines a discrete taylor-fourier transform problem.

Fields:
- `s::Vector{<:Real}`       | discrete signal [?]
- `t::StepRangeLen{<:Real}` | discrete time [s]
- `h::Vector{<:Int}`        | harmonic numbers [-]
- `D::Int`                  | degree of the highest derivative [-]
- `F::Real`                 | fundamental frequency [Hz]
- `K::Int`                  | o-spline degree [-]
- `N::Int`                  | number of samples per fundamental cycle [-]
"""
struct DTFTProblem <: AbstractDTFTProblem
    """"discrete signal"""
    s::Vector{<:Real}
    """discrete time"""
    t::StepRangeLen{<:Real}
    """harmonic numbers"""
    h::Vector{<:Int}
    """degree of the highest derivative"""
    D::Int
    """fundamental frequency"""
    F::Real
    """o-spline degree"""
    K::Int
    """number of samples per fundamental cycle"""
    N::Int

    # base constructor
    DTFTProblem(s::Vector{<:Real}, t::StepRangeLen{<:Real}, h::Vector{<:Int}, 
                D::Int, F::Real, K::Int, N::Int) =
        new(s, t, h, D, F, K, N)
end

# problem builder
function build_problem(s, t, h, D, F, K)
    # checks
    length(s) == length(t)          || error("the length of the discrete signal and time are inconsistent")
    all(diff(t) .≈ (t[2] - t[1]))   || error("the discrete time is not equidistance")

    # reduce the discrete time to a range `rT`
    # NB: the rounding function is added to avoid numerical issues
    ΔT  = round(t[2] - t[1], digits=10)
    rT  = round(t[1], digits=10):ΔT:round(t[end], digits=10)
    
    # determine the fundamental period [s]
    T   = 1 / F
    # find the number of samples in a fundamental cycle [-]
    # NB: the rounding function is added to avoid numerical issues
    N   = findfirst(x -> x == round(t[1] + T, digits=10), t) - 1

    return DTFTProblem(s, rT, h, D, F, K, N)
end
    
# solution types
"""
    TFT.AbstractDTFTSolution

Base for types which define discrete taylor-fourier transform solutions.
"""
abstract type AbstractDTFTSolution end

"""
    TFT.DTFTSolution

Defines a discrete taylor-fourier transform solution

Fields:
- `X::Dict{<:Int,Matrix{<:Complex}}`    | dictionary of dynamic phasors ξₕ⁽ᴰ⁾(t)
- `prob::TFT.AbstractDTFTProblem`       | DTFT problem struct
"""
struct DTFTSolution <: AbstractDTFTSolution 
    """dictionary of dynamic phasors ξₕ⁽ᴰ⁾"""
    X::Dict
    """DTFT problem struct"""
    prob::AbstractDTFTProblem
end