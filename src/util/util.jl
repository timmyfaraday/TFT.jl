################################################################################
#  Copyright 2022, Tom Van Acker (BASF Antwerp)                                #
################################################################################
# TFT.jl                                                                       #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TFT.jl                                    #
################################################################################

### TODO: frequency, rocof, error, checks

# checks
function check_sol(sol, D)
    sol.prob.D >= D || error("")
end

function check_sol(sol, D, H)
    sol.prob.D >= D || error("")
    H in sol.prob.h || error("")
end

# amplitude
"""
    TFT.a(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Shorthand function to obtain the amplitude of the Dth-degree derivative of the
Hth-harmonic phasor, dispatching to `amplitude(sol, D, H)`.
"""
a(sol::AbstractDTFTSolution, D::Int=0, H::Int=1) = amplitude(sol, D, H)

"""
    TFT.amplitude(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Function to obtain the amplitude of the Dth-degree derivative of the 
Hth-harmonic phasor.

See: [Assessing Synchrophasor Estimates of an Event Captured by a Phasor 
Measurement Unit, pg. 3112](https://ieeexplore.ieee.org/document/9239915)

Input:
- `sol::AbstractDFTFSolution`   | DTFT solution struct [-]
- `D::Int`                      | degree of the derivative [-], default=0
- `H::Int`                      | harmonic number [-], default=1

Output:
- `a::Vector{<:Real}`           | amplitude [?]
"""
function amplitude(sol::AbstractDTFTSolution, D::Int=0, H::Int=1)
    check_sol(sol, D, H)

    D == 0 && return 2.0 .* abs.(sol.X[H][:,D+1])
    D == 1 && return 2.0 .* real.(sol.X[H][:,D+1] .* exp.(-im .* φ(sol,0,H)))
    D == 2 && return 2.0 .* real.(sol.X[H][:,D+1] .* exp.(-im .* φ(sol,0,H))) .+ 
                        a(sol,0,H) .* φ(sol,1,H).^2
    
    return nothing
end

# angle
"""
    TFT.φ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Shorthand function to obtain the angle of the Dth-degree derivative of the
Hth-harmonic phasor, dispatching to `angle(sol, D, H)`.
"""
φ(sol::AbstractDTFTSolution, D::Int=0, H::Int=1) = angle(sol, D, H)

"""
    TFT.angle(TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Function to obtain the angle of the Dth-degree derivative of the 
Hth-harmonic phasor.
    
See: [Assessing Synchrophasor Estimates of an Event Captured by a Phasor 
Measurement Unit, pg. 3112](https://ieeexplore.ieee.org/document/9239915)
    
Input:
- `sol::AbstractDFTFSolution`   | DTFT solution struct [-]
- `D::Int`                      | degree of the derivative [-], default=0
- `H::Int`                      | harmonic number [-], default=1
    
Output:
- `φ::Vector{<:Real}`           | angle [(rad)]
"""
function angle(sol::AbstractDTFTSolution, D::Int=0, H::Int=1)
    check_sol(sol, D, H)

    # antirotated cte
    rS  = 0:(length(sol.prob.s)-1)
    ar  = exp.((-2.0 * pi * im * H / sol.prob.N) .* rS)

    D == 0 && return Base.angle.(sol.X[H][:,D+1] .* ar) 
    D == 1 && return imag.(sol.X[H][:,D+1] .* exp.(-im .* φh(sol,0,H))) ./
                        a(sol,0,H)
    D == 2 && return 2.0 .* imag.(sol.X[H][:,D+1] .* exp.(-im .* φh(sol,0,H))) ./
                        a(sol,0,H) .- a(sol,1,H) .* φh(sol,1,H)
    
    return nothing
end

# angle - alternative
"""
    TFT.φ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Shorthand function to obtain the angle of the Dth-degree derivative of the
Hth-harmonic phasor, dispatching to `angle(sol, D, H)`.
"""
φh(sol::AbstractDTFTSolution, D::Int=0, H::Int=1) = angleh(sol, D, H)

"""
    TFT.angle(TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Function to obtain the angle of the Dth-degree derivative of the 
Hth-harmonic phasor.
    
See: [Assessing Synchrophasor Estimates of an Event Captured by a Phasor 
Measurement Unit, pg. 3112](https://ieeexplore.ieee.org/document/9239915)
    
Input:
- `sol::AbstractDFTFSolution`   | DTFT solution struct [-]
- `D::Int`                      | degree of the derivative [-], default=0
- `H::Int`                      | harmonic number [-], default=1
    
Output:
- `φ::Vector{<:Real}`           | angle [(rad)]
"""
function angleh(sol::AbstractDTFTSolution, D::Int=0, H::Int=1)
    check_sol(sol, D, H)

    D == 0 && return Base.angle.(sol.X[H][:,D+1])
    D == 1 && return 2.0 .* imag.(sol.X[H][:,D+1] .* exp.(-im .* φh(sol,0,H))) ./
                        a(sol,0,H)
    D == 2 && return 2.0 .* imag.(sol.X[H][:,D+1] .* exp.(-im .* φh(sol,0,H))) ./
                        a(sol,0,H) .- a(sol,1,H) .* φh(sol,1,H)
    
    return nothing
end

# phasor
"""
    TFT.phasor(TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)
˚
Function to obtain the Dth-degree derivative of the Hth-harmonic phasor.

See: [Assessing Synchrophasor Estimates of an Event Captured by a Phasor 
Measurement Unit, pg. 3112](https://ieeexplore.ieee.org/document/9239915)
    
Input:
- `sol::AbstractDFTFSolution`   | DTFT solution struct [-]
- `D::Int`                      | degree of the derivative [-], default=0
- `H::Int`                      | harmonic number [-], default=1
    
Output:
- `phasor::Vector{<:Complex}`   | phasor [?]
"""
function phasor(sol::AbstractDTFTSolution, D::Int=0, H::Int=1)
    check_sol(sol, D, H)

    ω = 2.0 * pi * sol.prob.F * H

    D == 0 && return a(sol,0,H) .* exp.(im .* φ(sol,0,H))
    D == 1 && return (a(sol,1,H) .+ im .* a(sol,0,H) .* φ(sol,1,H)) .* 
                        exp.(im .* φ(sol,0,H))

    return nothing
end

# signal
"""
    TFT.signal(TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Function to obtain the Dth-degree derivative of the Hth-harmonic signal.

See: [Assessing Synchrophasor Estimates of an Event Captured by a Phasor 
Measurement Unit, pg. 3112](https://ieeexplore.ieee.org/document/9239915)
        
Input:
- `sol::AbstractDFTFSolution`   | DTFT solution struct [-]
- `D::Int`                      | degree of the derivative [-], default=0
- `H::Int`                      | harmonic number [-], default=1
        
Output:
- `signal::Vector{<:Real}`      | signal [?]
"""
function signal(sol::AbstractDTFTSolution, D::Int=0, H::Int=1)
    check_sol(sol, D, H)

    ω = 2.0 * pi * sol.prob.F * H

    D == 0 && return real(phasor(sol,0,H) .* exp.(im .* ω .* sol.prob.t))
    D == 1 && return real((phasor(sol,1,H) .+ im .* ω .* phasor(sol,0,H)) .*
                        exp.(im .* ω .* sol.prob.t))

    return nothing
end

# overall signal
"""
    TFT.overall_signal(TFT.AbstractDTFTSolution, D::Int=0)

Function to obtain the Dth-degree derivative of the overall signal.

See: [Assessing Synchrophasor Estimates of an Event Captured by a Phasor 
Measurement Unit, pg. 3112](https://ieeexplore.ieee.org/document/9239915)
        
Input:
- `sol::AbstractDFTFSolution`   | DTFT solution struct [-]
- `D::Int`                      | degree of the derivative [-], default=0
        
Output:
- `signal::Vector{<:Real}`      | signal [?]
"""
overall_signal(sol::AbstractDTFTSolution, D::Int=0) =
    sum(signal(sol,D,nh) for nh in sol.prob.h)

# error
"""
    TFT.error(sol::TFT.AbstractDTFTSolution)

Function to obtain the error between the input and computed signal.

Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]

Output:
- `e::Vector{<:Real}`           | error
"""
error(sol::TFT.AbstractDTFTSolution) = 
    sol.prob.s .- overall_signal(sol,0)