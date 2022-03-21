################################################################################
#  Copyright 2022, Tom Van Acker (BASF Antwerp)                                #
################################################################################
# TFT.jl                                                                       #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TFT.jl                                    #
################################################################################

# checks
function check_sol(sol, D)
    sol.prob.D >= D || Base.error("the required Dth-degree derivative is unavailable")
end

function check_sol(sol, D, H)
    sol.prob.D >= D || Base.error("the required Dth-degree derivative is unavailable")
    H in sol.prob.h || Base.error("the required Hth-harmonic phasor is unavailable")
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

    aₕ⁽⁰⁾(t) = |pₕ⁽⁰⁾(t)|
    aₕ⁽¹⁾(t) = ℜ[pₕ⁽¹⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] 
    aₕ⁽²⁾(t) = ℜ[pₕ⁽²⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] + aₕ⁽⁰⁾(t) ⋅ φₕ⁽¹⁾(t)^2

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

    D == 0 && return abs.(ξ(sol,D,H))
    D == 1 && return real.(ξ(sol,D,H) .* exp.(-im .* ϕ(sol,H)))
    D == 2 && return real.(ξ(sol,D,H) .* exp.(-im .* ϕ(sol,H))) .+ 
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

    φₕ⁽⁰⁾(t) = ∠[pₕ⁽⁰⁾(t) ⋅ exp(-2.0 ⋅ im ⋅ pi ⋅ H / N) ⋅ 0:length(s(t))-1]
    φₕ⁽¹⁾(t) = ℑ[pₕ⁽¹⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] / aₕ⁽⁰⁾(t)
    φₕ⁽²⁾(t) = (ℑ[pₕ⁽²⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] + 2 ⋅ aₕ⁽¹⁾(t) ⋅ φₕ⁽¹⁾(t)) / aₕ⁽⁰⁾(t)
    
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

    D == 0 && return Base.angle.(ξ(sol,D,H) .* ar) 
    D == 1 && return imag.(ξ(sol,D,H) .* exp.(-im .* ϕ(sol,H))) ./ a(sol,0,H)
    D == 2 && return (imag.(ξ(sol,D,H) .* exp.(-im .* ϕ(sol,H))) .- 
                        2.0 .* a(sol,1,H) .* φ(sol,1,H)
                     ) ./ a(sol,0,H)
    
    return nothing
end

# rotating angle
"""
    TFT.ϕ(sol::TFT.AbstractDTFTSolution, H::Int=1)

Shorthand function to obtain the alternative angle of the Dth-degree derivative 
of the Hth-harmonic phasor, dispatching to `angle(sol, D, H)`.
"""
ϕ(sol::AbstractDTFTSolution, H::Int=1) = rotating_angle(sol, H)

"""
    TFT.rotating_angle(TFT.AbstractDTFTSolution, H::Int=1)

Function to obtain the alternative angle of the Dth-degree derivative of the 
Hth-harmonic phasor.
    
    ϕₕ⁽⁰⁾(t) = ∠[pₕ⁽⁰⁾(t)]

See: [Assessing Synchrophasor Estimates of an Event Captured by a Phasor 
Measurement Unit, pg. 3112](https://ieeexplore.ieee.org/document/9239915)
    
Input:
- `sol::AbstractDFTFSolution`   | DTFT solution struct [-]
- `H::Int`                      | harmonic number [-], default=1
    
Output:
- `ϕ::Vector{<:Real}`           | rotating angle [(rad)]
"""
rotating_angle(sol::AbstractDTFTSolution, H::Int=1) =
    Base.angle.(ξ(sol,0,H))

# frequency
"""
    TFT.frequency(TFT.AbstractDTFTSolution, H::Int=1)

Function to obtain the Hth-harmonic frequency.

    fₕ(t) = Fₕ + φₕ⁽¹⁾(t) / (2 ⋅ π)

See: [new paper](xxx)

Input:
- `sol::AbstractDFTFSolution`   | DTFT solution struct [-]
- `H::Int`                      | harmonic number [-], default=1
    
Output:
- `freq::Vector{<:Real}`        | frequency [Hz/s]
"""
frequency(sol::AbstractDTFTSolution, H::Int=1) = 
    (H * sol.prob.F) .+ φ(sol,1,H) ./ (2 * pi)

# rocof
"""
    TFT.rocof(TFT.AbstractDTFTSolution, H::Int=1)

Function to obtain the Hth-harmonic rate-of-change-of-frequency.

See: [new paper](xxx)

    rₕ(t) = φₕ⁽²⁾(t) / (2 ⋅ π)^2

Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]
- `H::Int`                      | harmonic number [-], default=1

Output:
-`rocof::Vector{<:Real}`        | rate-of-change-of-frequency [Hz/s²]
"""
rocof(sol::AbstractDTFTSolution, H::Int=1) = φ(sol,2,H) / (2 * pi)^2

# complex_envelope
"""
    TFT.ξ(TFT.AbstractDFTFSolution, D::Int=0, H::Int=1)
"""
ξ(sol::AbstractDTFTSolution, D::Int=0, H::Int=1) = 2.0 .* sol.X[H][:,D+1]

# phasor
"""
    TFT.phasor(TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Function to obtain the Dth-degree derivative of the Hth-harmonic phasor.

    pₕ⁽ᴰ⁾(t) = 2 .* ξₕ⁽ᴰ⁾(t)

See: [xxx](xxx)
    
Input:
- `sol::AbstractDFTFSolution`   | DTFT solution struct [-]
- `D::Int`                      | degree of the derivative [-], default=0
- `H::Int`                      | harmonic number [-], default=1
    
Output:
- `phasor::Vector{<:Complex}`   | phasor [?]
"""
function phasor(sol::AbstractDTFTSolution, D::Int=0, H::Int=1)
    check_sol(sol, D, H)

    D == 0 && return a(sol,0,H) .* exp.(im .* φ(sol,0,H))
    D == 1 && return (a(sol,1,H) .+ im .* φ(sol,1,H) .* a(sol,0,H)
                     ) .* exp.(im .* φ(sol,0,H))
    D == 2 && return ((a(sol,2,H) .- φ(sol,1,H).^2 .* a(sol,0,H)) .+
                      im .* (2.0 .* φ(sol,1,H) .* a(sol,1,H) .+ φ(sol,2,H) .* a(sol,0,H))
                     ) .* exp.(im .* φ(sol,0,H))

    return nothing
end

# signal
"""
    TFT.signal(TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Function to obtain the Dth-degree derivative of the Hth-harmonic signal.

See: [Assessing Synchrophasor Estimates of an Event Captured by a Phasor 
Measurement Unit, pg. 3112](https://ieeexplore.ieee.org/document/9239915)
        
    sₕ⁽⁰⁾(t) = ℜ[pₕ⁽⁰⁾(t) ⋅ exp(im ⋅ ω ⋅ t)]
    sₕ⁽¹⁾(t) = ℜ[(pₕ⁽¹⁾(t) + im ⋅ ω ⋅ pₕ⁽⁰⁾(t)) ⋅ exp(im ⋅ ω ⋅ t)]
    sₕ⁽²⁾(t) = ℜ[(pₕ⁽²⁾(t) - ω² ⋅ pₕ⁽⁰⁾(t) + im ⋅ 2 ⋅ ω ⋅ pₕ⁽¹⁾(t)) ⋅ exp(im ⋅ ω ⋅ t)]

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
    D == 1 && return real(exp.(im .* ω .* sol.prob.t) .* 
                            (phasor(sol,1,H) .+ 
                             im .* ω .* phasor(sol,0,H))
                         )
    D == 2 && return real(exp.(im .* ω .* sol.prob.t) .* 
                            (phasor(sol,2,H) .- 
                             ω^2 .* phasor(sol,0,H) .+ 
                             im .* 2.0 .* ω .* phasor(sol,1,H)) 
                         )

    return nothing
end

# overall signal
"""
    TFT.overall_signal(TFT.AbstractDTFTSolution, D::Int=0)

Function to obtain the Dth-degree derivative of the overall signal.

    s⁽ᴰ⁾(t) = ∑ₕ sₕ⁽ᴰ⁾(t)

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

    e(t) = s(t) - s⁽⁰⁾(t)

Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]

Output:
- `e::Vector{<:Real}`           | error
"""
error(sol::TFT.AbstractDTFTSolution) = sol.prob.s .- overall_signal(sol,0)