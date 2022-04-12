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

# frequency and angular frequency
F(sol,H) = H * sol.prob.F
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

    aₕ⁽⁰⁾(t) = |2 ⋅ ξₕ⁽⁰⁾(t)| ∈ 𝐑⁺
    aₕ⁽¹⁾(t) = ℜ[2 ⋅ ξₕ⁽¹⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] ∈ 𝐑⁺
    aₕ⁽²⁾(t) = ℜ[2 ⋅ ξₕ⁽²⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] + aₕ⁽⁰⁾(t) ⋅ ϕₕ⁽¹⁾(t)² ∈ 𝐑⁺

See: [Assessing Synchrophasor Estimates of an Event Captured by a Phasor 
Measurement Unit, pg. 3112](https://ieeexplore.ieee.org/document/9239915)

Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]
- `D::Int`                      | degree of the derivative [-], default=0
- `H::Int`                      | harmonic number [-], default=1

Output:
- `a::Vector{<:Real}`           | amplitude aₕ⁽ᴰ⁾(t) ∈ 𝐑⁺ [?]
"""
function amplitude(sol::AbstractDTFTSolution, D::Int=0, H::Int=1)
    D == 0 && return abs.(2.0 .* ξ(sol,0,H))
    D == 1 && return real.(2.0 .* ξ(sol,1,H) .* exp.(-im .* ϕ(sol,0,H)))
    D == 2 && return real.(2.0 .* ξ(sol,2,H) .* exp.(-im .* ϕ(sol,0,H))) .+
                        a(sol,0,H) .* ϕ(sol,1,H).^2
    
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
    
    ϕₕ⁽⁰⁾(t) = ∠[pₕ⁽⁰⁾(t)] ∈ [-π,π]
    ϕₕ⁽¹⁾(t) = ℑ[2 ⋅ ξₕ⁽¹⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] / aₕ⁽⁰⁾(t) ∈ [-π,π]
    ϕₕ⁽²⁾(t) = {ℑ[2 ⋅ ξₕ⁽²⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] - 2 ⋅ aₕ⁽¹⁾(t) ⋅ ϕₕ⁽¹⁾(t)} / aₕ⁽⁰⁾(t) ∈ [-π,π]

See: [Assessing Synchrophasor Estimates of an Event Captured by a Phasor 
Measurement Unit, pg. 3112](https://ieeexplore.ieee.org/document/9239915)
    
Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]
- `D::Int`                      | degree of the derivative [-], default=0
- `H::Int`                      | harmonic number [-], default=1
    
Output:
- `ϕ::Vector{<:Real}`           | phase ϕₕ⁽ᴰ⁾(t) ∈ [-π,π] [(rad)]
"""
function phase(sol::AbstractDTFTSolution, D::Int=0, H::Int=1)
    D == 0 && return Base.angle.(ξ(sol,0,H))
    D == 1 && return imag.(2.0 .* ξ(sol,1,H) .* exp.(-im .* ϕ(sol,0,H))) ./
                        a(sol,0,H)
    D == 2 && return (imag.(2.0 .* ξ(sol,2,H) .* exp.(-im .* ϕ(sol,0,H))) .-
                        2.0 .* a(sol,1,H) .* ϕ(sol,1,H)) ./ a(sol,0,H)

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

    φₕ⁽⁰⁾(t) = ∠[ψₕ⁽⁰⁾(t)] ∈ [-π,π]
    φₕ⁽¹⁾(t) = ℑ[2 ⋅ ψₕ⁽¹⁾(t) ⋅ exp(-im ⋅ φₕ⁽⁰⁾(t))] / aₕ⁽⁰⁾(t) ∈ [-π,π]
    φₕ⁽²⁾(t) = {ℑ[2 ⋅ ψₕ⁽²⁾(t) ⋅ exp(-im ⋅ φₕ⁽⁰⁾(t))] - 2 ⋅ aₕ⁽¹⁾(t) ⋅ φₕ⁽¹⁾(t)} / aₕ⁽⁰⁾(t) ∈ [-π,π]
    
See: [Assessing Synchrophasor Estimates of an Event Captured by a Phasor 
Measurement Unit, pg. 3112](https://ieeexplore.ieee.org/document/9239915)
    
Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]
- `D::Int`                      | degree of the derivative [-], default=0
- `H::Int`                      | harmonic number [-], default=1
    
Output:
- `φ::Vector{<:Real}`           | anti-rotating phase φₕ⁽ᴰ⁾(t) ∈ [-π,π] [(rad)]
"""
function ar_phase(sol::AbstractDTFTSolution, D::Int=0, H::Int=1)
    D == 0 && return Base.angle.(ψ(sol,0,H)) 
    D == 1 && return imag.(2.0 .* ψ(sol,1,H) .* exp.(-im .* φ(sol,0,H))) ./ 
                        a(sol,0,H)
    D == 2 && return (imag.(2.0 .* ψ(sol,2,H) .* exp.(-im .* φ(sol,0,H))) .-
                        2.0 .* a(sol,1,H) .* φ(sol,1,H)) ./ a(sol,0,H)
    
    return nothing
end

# frequency
"""
    TFT.f(sol::TFT.AbstractDTFTSolution, H::Int=1)

Shorthand function to obtain the frequency of the zeroth-degree derivative of 
the Hth-harmonic phasor, dispatching to `frequency(sol,H)`.
"""
f(sol::AbstractDTFTSolution, H::Int=1) = frequency(sol,H)

"""
    TFT.frequency(sol::TFT.AbstractDTFTSolution, H::Int=1)

Function to obtain the frequency of the zeroth-degree derivative of the 
Hth-harmonic phasor.

    fₕ(t) = Fₕ + ϕₕ⁽¹⁾(t) / (2 π) ∈ 𝐑⁺

See: [Fast Taylor-Fourier Transform for Monitoring Modern Power Grids with 
Real-Time Dynamic Harmonic Estimation](tbp)

Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]
- `H::Int`                      | harmonic number [-], default=1
    
Output:
- `f::Vector{<:Real}`           | frequency fₕ(t) ∈ 𝐑⁺ [Hz]
"""
frequency(sol::AbstractDTFTSolution, H::Int=1) = 
    F(sol,H) .+ ϕ(sol,1,H) ./ (2 * pi)

# rocof
"""
    TFT.r(sol::TFT.AbstractDTFTSolution, H::Int=1)

Shorthand function to obtain the rate-of-change-of-frequency of the 
zeroth-degree derivative of the Hth-harmonic phasor, dispatching to 
`rocof(sol,H)`.
"""
r(sol::AbstractDTFTSolution, H::Int=1) = rocof(sol,H)

"""
    TFT.rocof(sol::TFT.AbstractDTFTSolution, H::Int=1)

Function to obtain the rate-of-change-of-frequency of the zeroth-degree 
derivative of the Hth-harmonic phasor.

    rₕ(t) = ϕₕ⁽²⁾(t) / (2 π)² ∈ 𝐑

See: [Fast Taylor-Fourier Transform for Monitoring Modern Power Grids with 
Real-Time Dynamic Harmonic Estimation](tbp)

Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]
- `H::Int`                      | harmonic number [-], default=1
    
Output:
- `r::Vector{<:Real}`           | rocof rₕ(t) ∈ 𝐑 [Hz/s]
"""
rocof(sol::AbstractDTFTSolution, H::Int=1) = 
    ϕ(sol,2,H) ./ (2 * pi)^2

# dynamic phasor
"""
    TFT.ξ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Shorthand function to obtain the Dth-degree derivative of the Hth-harmonic 
dynamic phasor, dispatching to `phasor(sol,D,H)`.
"""
ξ(sol::AbstractDTFTSolution, D::Int=0, H::Int=1) = phasor(sol,D,H)

"""
    TFT.phasor(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Function to obtain the D-th-degree derivative of the Hth-harmonic dynamic 
phasor.

Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]
- `D::Int`                      | degree of the derivative [-], default=0
- `H::Int`                      | harmonic number [-], default=1

Output:
- `ξ::Vector{<:Complex}`        | dynamic phasor ξₕ⁽ᴰ⁾(t) ∈ 𝐂 [?]
"""
function phasor(sol::AbstractDTFTSolution, D::Int=0, H::Int=1) 
    check_sol(sol, D, H)
    
    return sol.X[H][:,D+1]
end 

# anti-rotating dynamic phasor
"""
    TFT.ψ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Shorthand function to obtain the Dth-degree derivative of the Hth-harmonic anti-
rotating dynamic phasor, dispatching to `ar_phasor(sol,D,H)`.
"""
ψ(sol::AbstractDTFTSolution, D::Int=0, H::Int=1) = 
    ar_phasor(sol,D,H)

"""
    TFT.ar_phasor(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)

Function to obtain the D-th-degree derivative of the Hth-harmonic anti-rotating 
dynamic phasor.

    ψₕ⁽ᴰ⁾(t) = ξₕ⁽ᴰ⁾(t) exp(-im ωₕ t) ∈ 𝐂

Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]
- `D::Int`                      | degree of the derivative [-], default=0
- `H::Int`                      | harmonic number [-], default=1

Output:
- `ψ::Vector{<:Complex}`        | anti-rotating dynamic phasor ψₕ⁽ᴰ⁾(t) ∈ 𝐂 [?]
"""
ar_phasor(sol::AbstractDTFTSolution, D::Int=0, H::Int=1) =
    ξ(sol,D,H) .* exp.(-im .* ω(sol,H) .* sol.prob.t)

# signal
"""
    TFT.signal(sol::TFT.AbstractDTFTSolution)

Function to obtain the overall signal.

    s(t) = ∑ₕ ξₕ⁽⁰⁾(t) + conj(ξₕ⁽⁰⁾(t)) ∈ 𝐑

Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]
            
Output:
- `s::Vector{<:Real}`           | signal s(t) ∈ 𝐑 [?]
"""
signal(sol::AbstractDTFTSolution) = 
    sum(real(ξ(sol,0,nh) + conj(ξ(sol,0,nh))) for nh in sol.prob.h)

"""
    TFT.signal(sol::TFT.AbstractDTFTSolution, H::Int)

Function to obtain the Hth-harmonic signal. 

    sₕ(t) = ξₕ⁽⁰⁾(t) + conj(ξₕ⁽⁰⁾(t)) ∈ 𝐑

Input:
- `sol::AbstractDTFTSolution`   | DTFT solution struct [-]
- `H::Int`                      | harmonic number [-]
        
Output:
- `s::Vector{<:Real}`           | signal sₕ(t) ∈ 𝐑 [?]
"""
signal(sol::AbstractDTFTSolution, H::Int) = real(ξ(sol,0,H) + conj(ξ(sol,0,H)))

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
