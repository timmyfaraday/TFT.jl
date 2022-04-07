################################################################################
#  Copyright 2022, Tom Van Acker (BASF Antwerp)                                #
################################################################################
# TFT.jl                                                                       #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TFT.jl                                    #
################################################################################

# checks
function check_sol(sol, D, H)
    sol.prob.D >= D || Base.error("the required Dth-degree derivative is unavailable")
    H in sol.prob.h || Base.error("the required Hth-harmonic phasor is unavailable")
end

# angular frequency
ω(sol,H) = 2.0 * pi * H * sol.prob.F 

# amplitude
"""
    TFT.a(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Shorthand function to obtain the amplitude of the Dth-degree derivative of the
Hth-harmonic phasor, dispatching to `amplitude(sol,D,H)`.
"""
a(sol::AbstractDTFTSolution, D::Int=0, H::Int=1) = amplitude(sol,D,H)

"""
    TFT.amplitude(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Function to obtain the amplitude of the Dth-degree derivative of the 
Hth-harmonic phasor.

    aₕ⁽⁰⁾(t) = 2 ⋅ |ξₕ⁽⁰⁾(t)|
    # aₕ⁽¹⁾(t) = 2 ⋅ ℜ[ξₕ⁽¹⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] 
    # aₕ⁽²⁾(t) = 2 ⋅ ℜ[ξₕ⁽²⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] + aₕ⁽⁰⁾(t) ⋅ φₕ⁽¹⁾(t)^2

See: [Assessing Synchrophasor Estimates of an Event Captured by a Phasor 
Measurement Unit, pg. 3112](https://ieeexplore.ieee.org/document/9239915)

Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]
- `D::Int`                      | degree of the derivative [-], default=0
- `H::Int`                      | harmonic number [-], default=1

Output:
- `a::Vector{<:Real}`           | amplitude aₕ⁽ᴰ⁾(t) [?]
"""
function amplitude(sol::AbstractDTFTSolution, D::Int=0, H::Int=1)
    D == 0 && return 2.0 .* abs.(ξ(sol,D,H))
    # D == 1 && return 2.0 .* real.(ξ(sol,D,H) .* exp.(-im .* ϕ(sol,H)))
    # D == 2 && return 2.0 .* real.(ξ(sol,D,H) .* exp.(-im .* ϕ(sol,H))) .+ 
    #                    a(sol,0,H) .* φ(sol,1,H).^2
    
    return nothing
end

# phase
"""
    TFT.ϕ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Shorthand function to obtain the phase of the Dth-degree derivative 
of the Hth-harmonic phasor, dispatching to `phase(sol, D, H)`.
"""
ϕ(sol::AbstractDTFTSolution, D::Int=0, H::Int=1) = phase(sol, D, H)

"""
    TFT.phase(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Function to obtain the alternative angle of the Dth-degree derivative of the 
Hth-harmonic phasor.
    
    ϕₕ⁽⁰⁾(t) = ∠[pₕ⁽⁰⁾(t)]

See: [Assessing Synchrophasor Estimates of an Event Captured by a Phasor 
Measurement Unit, pg. 3112](https://ieeexplore.ieee.org/document/9239915)
    
Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]
- `H::Int`                      | harmonic number [-], default=1
    
Output:
- `ϕ::Vector{<:Real}`           | phase ϕₕ⁽ᴰ⁾(t) [(rad)]
"""
function phase(sol::AbstractDTFTSolution, D::Int=0, H::Int=1)
    D == 0 && return Base.angle.(ξ(sol,0,H))
    # D == 1 
    # D == 2

    return nothing
end

# anti-rotating phase
"""
    TFT.φ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Shorthand function to obtain the anti-rotating phase of the Dth-degree 
derivative of the Hth-harmonic phasor, dispatching to `angle(sol,D,H)`.
"""
φ(sol::AbstractDTFTSolution, D::Int=0, H::Int=1) = ar_phase(sol,D,H)

"""
    TFT.ar_phase(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Function to obtain the anti-rotating phase of the Dth-degree derivative of the 
Hth-harmonic phasor.

    φₕ⁽⁰⁾(t) = ∠[ξₕ⁽⁰⁾(t) ⋅ exp(-2.0 ⋅ im ⋅ pi ⋅ H / N) ⋅ 0:length(s(t))-1]
    # φₕ⁽¹⁾(t) = ℑ[pₕ⁽¹⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] / aₕ⁽⁰⁾(t)
    # φₕ⁽²⁾(t) = (ℑ[pₕ⁽²⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] + 2 ⋅ aₕ⁽¹⁾(t) ⋅ φₕ⁽¹⁾(t)) / aₕ⁽⁰⁾(t)
    
See: [Assessing Synchrophasor Estimates of an Event Captured by a Phasor 
Measurement Unit, pg. 3112](https://ieeexplore.ieee.org/document/9239915)
    
Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]
- `D::Int`                      | degree of the derivative [-], default=0
- `H::Int`                      | harmonic number [-], default=1
    
Output:
- `φ::Vector{<:Real}`           | anti-rotating phase φₕ⁽ᴰ⁾(t) [(rad)]
"""
function ar_phase(sol::AbstractDTFTSolution, D::Int=0, H::Int=1)
    # anti-rotating cte
    rS  = 0:(length(sol.prob.s)-1)
    ar  = exp.((-2.0 * pi * im * H / sol.prob.N) .* rS)

    D == 0 && return Base.angle.(ξ(sol,D,H) .* ar) 
    # D == 1 && return imag.(ξ(sol,D,H) .* exp.(-im .* ϕ(sol,H))) ./ a(sol,0,H)
    # D == 2 && return (imag.(ξ(sol,D,H) .* exp.(-im .* ϕ(sol,H))) .- 
    #                     2.0 .* a(sol,1,H) .* φ(sol,1,H)
    #                  ) ./ a(sol,0,H)
    
    return nothing
end

# # frequency
# """
#     TFT.frequency(TFT.AbstractDTFTSolution, H::Int=1)

# Function to obtain the Hth-harmonic frequency.

#     fₕ(t) = Fₕ + φₕ⁽¹⁾(t) / (2 ⋅ π)

# See: [new paper](xxx)

# Input:
# - `sol::AbstractDTFTSolution`   | DTFT solution struct [-]
# - `H::Int`                      | harmonic number [-], default=1
    
# Output:
# - `freq::Vector{<:Real}`        | frequency f(t) [Hz/s]
# """
# frequency(sol::AbstractDTFTSolution, H::Int=1) = 
#     (H * sol.prob.F) .+ φ(sol,1,H) ./ (2 * pi)

# # rocof
# """
#     TFT.rocof(TFT.AbstractDTFTSolution, H::Int=1)

# Function to obtain the Hth-harmonic rate-of-change-of-frequency.

# See: [new paper](xxx)

#     rₕ(t) = φₕ⁽²⁾(t) / (2 ⋅ π)^2

# Input:
# - `sol::AbstractDTFTSolution`   | DTFT solution struct [-]
# - `H::Int`                      | harmonic number [-], default=1

# Output:
# -`rocof::Vector{<:Real}`        | rate-of-change-of-frequency r(t) [Hz/s²]
# """
# rocof(sol::AbstractDTFTSolution, H::Int=1) = φ(sol,2,H) / (2 * pi)^2

# dynamic phasor
"""
    TFT.ξ(TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Shorthand function to obtain the Dth-degree derivative of the Hth-harmonic 
dynamic phasor, dispatching to `phasor(sol,D,H)`.
"""
ξ(sol::AbstractDTFTSolution, D::Int=0, H::Int=1) = phasor(sol,D,H)

"""
    TFT.phasor(TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Function to obtain the D-th-degree derivative of the Hth-harmonic dynamic 
phasor.

Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]
- `D::Int`                      | degree of the derivative [-], default=0
- `H::Int`                      | harmonic number [-], default=1

Output:
- phasor::Vector{<:Complex}     | dynamic phasor ξₕ⁽ᴰ⁾(t) [?]
"""
function phasor(sol::AbstractDTFTSolution, D::Int=0, H::Int=1) 
    check_sol(sol, D, H)
    
    return sol.X[H][:,D+1]
end 

# anti-rotating dynamic phasor
"""
    TFT.ψ(TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Shorthand function to obtain the Dth-degree derivative of the Hth-harmonic anti-
rotating dynamic phasor, dispatching to `ar_phasor(sol,D,H)`.
"""
ψ(sol::AbstractDTFTSolution, D::Int=0, H::Int=1) = 
    ar_phasor(sol,D,H)

"""
    TFT.ar_phasor(TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Function to obtain the D-th-degree derivative of the Hth-harmonic anti-rotating 
dynamic phasor.

    ψₕ⁽ᴰ⁾(t) = ξₕ⁽ᴰ⁾(t) exp(-im ωₕ t)

Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]
- `D::Int`                      | degree of the derivative [-], default=0
- `H::Int`                      | harmonic number [-], default=1

Output:
- ar_phasor::Vector{<:Complex}  | anti-rotating dynamic phasor ψₕ⁽ᴰ⁾(t) [?]
"""
ar_phasor(sol::AbstractDTFTSolution, D::Int=0, H::Int=1) =
    ξ(sol,D,H) .* exp.(-im .* ω(sol,H) .* sol.prob.t)

# signal
"""
    TFT.signal(TFT.AbstractDTFTSolution)

Function to obtain the overall signal.

    s(t) = ∑ₕ ξₕ⁽⁰⁾(t) + conj(ξₕ⁽⁰⁾(t))

Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]
            
Output:
- `signal::Vector{<:Real}`      | signal s(t) [?]
"""
signal(sol::AbstractDTFTSolution) = 
    sum(real(ξ(sol,0,nh) + conj(ξ(sol,0,nh))) for nh in sol.prob.h)

"""
    TFT.signal(TFT.AbstractDTFTSolution, H::Int)

Function to obtain the Hth-harmonic signal. 

    sₕ(t) = ξₕ⁽⁰⁾(t) + conj(ξₕ⁽⁰⁾(t))

Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]
- `H::Int`                      | harmonic number [-]
        
Output:
- `signal::Vector{<:Real}`      | signal sₕ(t) [?]
"""
signal(sol::AbstractDTFTSolution, H::Int) =
    real(ξ(sol,0,H) + conj(ξ(sol,0,H)))

# error
"""
    TFT.error(sol::TFT.AbstractDTFTSolution)

Function to obtain the error between the input and computed signal.

Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]

Output:
- `e::Vector{<:Real}`           | error [?]
"""
error(sol::TFT.AbstractDTFTSolution) = sol.prob.s .- signal(sol)
