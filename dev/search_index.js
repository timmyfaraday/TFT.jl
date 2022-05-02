var documenterSearchIndex = {"docs":
[{"location":"guide/#TFT-coding-style-guide","page":"Developer's Guide","title":"TFT coding style guide","text":"","category":"section"},{"location":"guide/","page":"Developer's Guide","title":"Developer's Guide","text":"a constant is denoted by a capital letter,\na vector is denoted using a small letter\na matrix is denoted by a capital letter,\na range is denoted by r, followed by a capital letter,\na step is preceeded by a capital greek Δ","category":"page"},{"location":"dtft/#Discrete-Taylor-Fourier-Transform","page":"DTFT","title":"Discrete Taylor Fourier Transform","text":"","category":"section"},{"location":"dtft/#DTFT-Problem","page":"DTFT","title":"DTFT Problem","text":"","category":"section"},{"location":"dtft/","page":"DTFT","title":"DTFT","text":"TFT.AbstractDTFTProblem","category":"page"},{"location":"dtft/#TFT.AbstractDTFTProblem","page":"DTFT","title":"TFT.AbstractDTFTProblem","text":"TFT.AbstractDTFTProblem\n\nBase for types which define discrete taylor-fourier transform problems.\n\n\n\n\n\n","category":"type"},{"location":"dtft/","page":"DTFT","title":"DTFT","text":"TFT.DTFTProblem","category":"page"},{"location":"dtft/#TFT.DTFTProblem","page":"DTFT","title":"TFT.DTFTProblem","text":"TFT.DTFTProblem\n\nDefines a discrete taylor-fourier transform problem.\n\nFields:\n\ns::Vector{<:Real}       | discrete signal [?]\nt::StepRangeLen{<:Real} | discrete time [s]\nh::Vector{<:Int}        | harmonic numbers [-]\nD::Int                  | degree of the highest derivative [-]\nF::Real                 | fundamental frequency [Hz]\nK::Int                  | o-spline degree [-]\nN::Int                  | number of samples per fundamental cycle [-]\n\n\n\n\n\n","category":"type"},{"location":"dtft/#DTFT-Solution","page":"DTFT","title":"DTFT Solution","text":"","category":"section"},{"location":"dtft/","page":"DTFT","title":"DTFT","text":"TFT.AbstractDTFTSolution","category":"page"},{"location":"dtft/#TFT.AbstractDTFTSolution","page":"DTFT","title":"TFT.AbstractDTFTSolution","text":"TFT.AbstractDTFTSolution\n\nBase for types which define discrete taylor-fourier transform solutions.\n\n\n\n\n\n","category":"type"},{"location":"dtft/","page":"DTFT","title":"DTFT","text":"TFT.DTFTSolution","category":"page"},{"location":"dtft/#TFT.DTFTSolution","page":"DTFT","title":"TFT.DTFTSolution","text":"TFT.DTFTSolution\n\nDefines a discrete taylor-fourier transform solution\n\nFields:\n\nX::Dict{<:Int,Matrix{<:Complex}}    | dictionary of dynamic phasors ξₕ⁽ᴰ⁾(t)\nprob::TFT.AbstractDTFTProblem       | DTFT problem struct\n\n\n\n\n\n","category":"type"},{"location":"math/#Mathematical-Background-on-Taylor-Fourier-Transform","page":"Mathematical Background","title":"Mathematical Background on Taylor-Fourier Transform","text":"","category":"section"},{"location":"math/#Nomenclature","page":"Mathematical Background","title":"Nomenclature","text":"","category":"section"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"Symbol Description Domain\ns signal 𝐑\na amplitude 𝐑⁺\nϕ phase [(rad)] [-π,π]\nφ anti-rotating phase [(rad)] [-π,π]\nξ dynamic phasor 𝐂\nψ anti-rotating dynamic phasor 𝐂","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"Set Description Domain\nd ∈ D set of derivatives 𝐍\nh ∈ H set of harmonic number 𝐍","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"Parameter Description Domain\nF frequency [Hz] 𝐑⁺\nω angular frequency [(rad)/s] 𝐑⁺","category":"page"},{"location":"math/#Taylor-Fourier-Transform","page":"Mathematical Background","title":"Taylor-Fourier Transform","text":"","category":"section"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"The Taylor-Fourier Transform function tft() gives the up-to-Dth derivative of  the Hth-harmonic dynamic phasors ξₕ⁽ᵈ⁾(t) ∈ ℂ, ∀ d ∈ {0,..,D}, h ∈ H. ","category":"page"},{"location":"math/#Zeroth-Harmonic","page":"Mathematical Background","title":"Zeroth Harmonic","text":"","category":"section"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"The up-to-Dth derivative of the zeroth-harmonic dynamic phasor  ξ₀⁽ᵈ⁾(t) ∈ ℝ, ∀ d ∈ {0,..,D} is given by,","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"    xi^(d)_0(t)        = a^(d)_0(t)","category":"page"},{"location":"math/#Zeroth-Derivative","page":"Mathematical Background","title":"Zeroth Derivative","text":"","category":"section"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"The zeroth derivative of the hth-harmonic dynamic phasor ξₕ⁽⁰⁾(t) is given by,","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"beginaligned\n    xi^(0)_h(t)        = fraca^(0)_h(t)2  exp(j phi^(0)_h(t)) \n                            = fraca^(0)_h(t)2  exp(j varphi^(0)_h(t))  exp(j omega_h t) \n                            = psi^(0)_h(t)  exp(j omega_h t)\n\nendaligned","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"where ψₕ⁽⁰⁾(t) denotes the anti-rotating zeroth derivative of the hth-harmonic  dynamic phasor. ","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"The zeroth derivative of the amplitude aₕ⁽⁰⁾(t), phase ϕₕ⁽⁰⁾(t) and anti-rotating  phase φₕ⁽⁰⁾(t) are, respectively,","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"beginaligned\n    a^(0)_h(t)          = 2  xi^(0)_h(t) \n    phi^(0)_h(t)       = angle big( xi^(0)_h(t) big) \n    varphi^(0)_h(t)    = angle big( xi^(0)_h(t)  exp(-j omega_h t) big) = angle big( psi^(0)_h(t) big)\nendaligned","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"The zeroth derivative hth-harmonic signal sₕ⁽⁰⁾(t) and overall signal s⁽⁰⁾(t)  are, respectively,","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"beginaligned\n    s^(0)_h(t)          = xi^(0)_h(t) + textconjbig( xi^(0)_h(t) big) \n    s^(0)(t)              = sum_h in H xi^(0)_h(t) + textconjbig( xi^(0)_h(t) big)\nendaligned","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"where conj(ξₕ⁽⁰⁾(t)) denotes the complex conjugate of the zeroth derivative of the  hth-harmonic dynamic phasor.","category":"page"},{"location":"math/#First-Derivative","page":"Mathematical Background","title":"First Derivative","text":"","category":"section"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"The first derivate of the hth-harmonic dynamic phasor ξₕ⁽¹⁾(t) is given by,","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"beginaligned\n    xi_h^(1)(t)        =  fracmathrmdxi_h^(0)(t)mathrmdt \n                            =  frac12 big( a^(1)_h(t) +\n                                j  phi^(1)_h(t)  a^(0)_h(t)) big)\n                                 exp(j phi^(0)_h(t))\n\nendaligned","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"The first derivative of the hth-harmomic anti-rotating dynamic phasor ψₕ⁽¹⁾(t) is,","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"beginaligned\n    psi_h^(1)(t)       =  fracmathrmdpsi_h^(0)(t)mathrmdt \n                            =  xi_h^(1)(t)  exp(-j omega_h t)\nendaligned","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"The first derivative of the hth-harmonic amplitude aₕ⁽¹⁾(t), phase ϕₕ⁽¹⁾(t) and anti-rotating phase φₕ⁽¹⁾(t) are, respectively,","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"beginaligned\n    a_h^(1)(t)          = ℜ2  xi_h^(1)(t)  exp(-j phi_h^(0)(t)) \n    phi_h^(1)(t)       = fracℑ2  xi_h^(1)(t)  exp(-j phi_h^(0)(t))a^(0)_h(t) \n    varphi_h^(1)(t)    = fracℑ2  psi_h^(1)(t)  exp(-j varphi_h^(0)(t))a^(0)_h(t)\nendaligned","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"Note that implicitly this means that ϕₕ⁽¹⁾(t) = φₕ⁽¹⁾(t).","category":"page"},{"location":"math/#Second-Derivative","page":"Mathematical Background","title":"Second Derivative","text":"","category":"section"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"The second derivate of the hth-harmonic dynamic phasor ξₕ⁽²⁾(t) is given by,","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"beginaligned\n    xi_h^(2)(t)        =  fracmathrmdxi_h^(1)(t)mathrmdt \n                            =  frac12 big( \n                                a^(2)_h(t) - \n                                phi^(1)_h(t)^2  a^(0)_h(t) +\n                                j  2  phi^(1)_h(t)  a^(1)_h(t) +\n                                j  phi^(2)_h(t)  a^(0)_h(t)) big)\n                                 exp(j phi^(0)_h(t))\n\nendaligned","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"The second derivative of the hth-harmomic anti-rotating dynamic phasor ψₕ⁽²⁾(t) is,","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"beginaligned\n    psi_h^(2)(t)       =  fracmathrmdpsi_h^(1)(t)mathrmdt \n                            =  xi_h^(2)(t) exp(-j omega_h t)\nendaligned","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"The second derivative of the hth-harmonic amplitude aₕ⁽²⁾(t), phase ϕₕ⁽²⁾(t) and anti-rotating phase φₕ⁽²⁾(t) are, respectively,","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"beginaligned\n    a_h^(2)(t)          = ℜ2  xi_h^(2)(t)  exp(-j phi_h^(0)(t)) + a_h^(0)(t)  phi_h^(1)(t)^2 \n    phi_h^(1)(t)       = fracℑ2  xi_h^(2)(t)  exp(-j phi_h^(0)(t)) - 2  a_h^(1)(t)  phi_h^(1)(t)a^(0)_h(t) \n    varphi_h^(1)(t)    = fracℑ2  psi_h^(2)(t)  exp(-j varphi_h^(0)(t)) - 2  a_h^(1)(t)  varphi_h^(1)(t)a^(0)_h(t)\nendaligned","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"Note that implicitly this means that ϕₕ⁽²⁾(t) = φₕ⁽²⁾(t).","category":"page"},{"location":"#TFT.jl","page":"Home","title":"TFT.jl","text":"","category":"section"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"TFT.jl is a Julia package for Taylor-Fourier Transform.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The latest stable release of TFT.jl can be installed using the Julia package  manager:","category":"page"},{"location":"","page":"Home","title":"Home","text":"] add \"https://github.com/timmyfaraday/TFT.jl\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"In order to test whether the package works, run:","category":"page"},{"location":"","page":"Home","title":"Home","text":"] test TFT","category":"page"},{"location":"util/#Utilities","page":"Utilities","title":"Utilities","text":"","category":"section"},{"location":"util/","page":"Utilities","title":"Utilities","text":"A number of functions are made available to the user to retrieve specific components of the TFT solution","category":"page"},{"location":"util/#Amplitude","page":"Utilities","title":"Amplitude","text":"","category":"section"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.a(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.a","page":"Utilities","title":"TFT.a","text":"TFT.a(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nShorthand function to obtain the amplitude of the Dth-degree derivative of the Hth-harmonic phasor, dispatching to amplitude(sol,D,H).\n\n\n\n\n\n","category":"function"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.amplitude(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.amplitude","page":"Utilities","title":"TFT.amplitude","text":"TFT.amplitude(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nFunction to obtain the amplitude of the Dth-degree derivative of the  Hth-harmonic phasor.\n\n∀ h ∈ {0}:\n    a₀⁽ᴰ⁾(t) = ξ₀⁽ᴰ⁾(t)\n∀ h ∈ 𝓗/{0}:\n    aₕ⁽⁰⁾(t) = |2 ⋅ ξₕ⁽⁰⁾(t)| ∈ 𝐑⁺\n    aₕ⁽¹⁾(t) = ℜ[2 ⋅ ξₕ⁽¹⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] ∈ 𝐑\n    aₕ⁽²⁾(t) = ℜ[2 ⋅ ξₕ⁽²⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] + aₕ⁽⁰⁾(t) ⋅ ϕₕ⁽¹⁾(t)² ∈ 𝐑\n\nSee: Assessing Synchrophasor Estimates of an Event Captured by a Phasor  Measurement Unit, pg. 3112\n\nInput:\n\nsol::AbstractDTFTSolution   | DTFT solution struct [-]\nD::Int                      | degree of the derivative [-], default=0\nH::Int                      | harmonic number [-], default=1\n\nOutput:\n\na::Vector{<:Real}           | amplitude aₕ⁽ᴰ⁾(t) [?]\n\n\n\n\n\n","category":"function"},{"location":"util/","page":"Utilities","title":"Utilities","text":"note: Note\nThe amplitude aₕ⁽¹⁾(t) denotes the first derivate of the amplitude of the  zeroth derivative of the dynamic phasor ξₕ⁽⁰⁾(t), not the amplitude of the  first derivative of the dynamic phasor ξₕ⁽¹⁾(t).","category":"page"},{"location":"util/#Phase","page":"Utilities","title":"Phase","text":"","category":"section"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.ϕ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.ϕ","page":"Utilities","title":"TFT.ϕ","text":"TFT.ϕ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nShorthand function to obtain the phase of the Dth-degree derivative  of the Hth-harmonic phasor, dispatching to phase(sol, D, H).\n\n\n\n\n\n","category":"function"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.phase(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.phase","page":"Utilities","title":"TFT.phase","text":"TFT.phase(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nFunction to obtain the alternative angle of the Dth-degree derivative of the  Hth-harmonic phasor.\n\n∀ h in {0}:\n    ϕ₀⁽ᴰ⁾(t) = 0.0\n∀ h ∈ 𝓗/{0}:\n    ϕₕ⁽⁰⁾(t) = ∠[pₕ⁽⁰⁾(t)] ∈ [-π,π]\n    ϕₕ⁽¹⁾(t) = ℑ[2 ⋅ ξₕ⁽¹⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] / aₕ⁽⁰⁾(t) ∈ 𝐑\n    ϕₕ⁽²⁾(t) = {ℑ[2 ⋅ ξₕ⁽²⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] - 2 ⋅ aₕ⁽¹⁾(t) ⋅ ϕₕ⁽¹⁾(t)} / aₕ⁽⁰⁾(t) ∈ 𝐑\n\nSee: Assessing Synchrophasor Estimates of an Event Captured by a Phasor  Measurement Unit, pg. 3112\n\nInput:\n\nsol::AbstractDTFTSolution   | DTFT solution struct [-]\nD::Int                      | degree of the derivative [-], default=0\nH::Int                      | harmonic number [-], default=1\n\nOutput:\n\nϕ::Vector{<:Real}           | phase ϕₕ⁽ᴰ⁾(t) [(rad)]\n\n\n\n\n\n","category":"function"},{"location":"util/","page":"Utilities","title":"Utilities","text":"note: Note\nThe phase ϕₕ⁽¹⁾(t) denotes the first derivate of the phase of the  zeroth derivative of the dynamic phasor ξₕ⁽⁰⁾(t), not the phase of the first  derivative of the dynamic phasor ξₕ⁽¹⁾(t).","category":"page"},{"location":"util/#Anti-Rotating-Phase","page":"Utilities","title":"Anti-Rotating Phase","text":"","category":"section"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.φ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.φ","page":"Utilities","title":"TFT.φ","text":"TFT.φ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nShorthand function to obtain the anti-rotating phase of the Dth-degree  derivative of the Hth-harmonic phasor, dispatching to angle(sol,D,H).\n\n\n\n\n\n","category":"function"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.ar_phase(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.ar_phase","page":"Utilities","title":"TFT.ar_phase","text":"TFT.ar_phase(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nFunction to obtain the anti-rotating phase of the Dth-degree derivative of the  Hth-harmonic phasor.\n\n∀ h in {0}:\n    φ₀⁽ᴰ⁾(t) = 0.0\n∀ h ∈ 𝓗/{0}:\n    φₕ⁽⁰⁾(t) = ∠[ψₕ⁽⁰⁾(t)] ∈ [-π,π]\n    φₕ⁽¹⁾(t) = ℑ[2 ⋅ ψₕ⁽¹⁾(t) ⋅ exp(-im ⋅ φₕ⁽⁰⁾(t))] / aₕ⁽⁰⁾(t) ∈ 𝐑\n    φₕ⁽²⁾(t) = {ℑ[2 ⋅ ψₕ⁽²⁾(t) ⋅ exp(-im ⋅ φₕ⁽⁰⁾(t))] - 2 ⋅ aₕ⁽¹⁾(t) ⋅ φₕ⁽¹⁾(t)} / aₕ⁽⁰⁾(t) ∈ 𝐑\n\nSee: Assessing Synchrophasor Estimates of an Event Captured by a Phasor  Measurement Unit, pg. 3112\n\nInput:\n\nsol::AbstractDTFTSolution   | DTFT solution struct [-]\nD::Int                      | degree of the derivative [-], default=0\nH::Int                      | harmonic number [-], default=1\n\nOutput:\n\nφ::Vector{<:Real}           | anti-rotating phase φₕ⁽ᴰ⁾(t) [(rad)]\n\n\n\n\n\n","category":"function"},{"location":"util/","page":"Utilities","title":"Utilities","text":"note: Note\nThe anti-rotating phase φₕ⁽¹⁾(t) denotes the first derivate of the  anti-rotating phase of the zeroth derivative of the dynamic phasor  ξₕ⁽⁰⁾(t), not the anti-rotating phase of the first derivative of the  dynamic phasor ξₕ⁽¹⁾(t).","category":"page"},{"location":"util/#Frequency","page":"Utilities","title":"Frequency","text":"","category":"section"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.f(sol::TFT.AbstractDTFTSolution, H::Int=1)","category":"page"},{"location":"util/#TFT.f","page":"Utilities","title":"TFT.f","text":"TFT.f(sol::TFT.AbstractDTFTSolution, H::Int=1)\n\nShorthand function to obtain the frequency of the zeroth-degree derivative of  the Hth-harmonic phasor, dispatching to frequency(sol,H).\n\n\n\n\n\n","category":"function"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.frequency(sol::TFT.AbstractDTFTSolution, H::Int=1)","category":"page"},{"location":"util/#TFT.frequency","page":"Utilities","title":"TFT.frequency","text":"TFT.frequency(sol::TFT.AbstractDTFTSolution, H::Int=1)\n\nFunction to obtain the frequency of the zeroth-degree derivative of the  Hth-harmonic phasor.\n\nfₕ(t) = Fₕ + ϕₕ⁽¹⁾(t) / (2 π) ∈ 𝐑⁺\n\nSee: Fast Taylor-Fourier Transform for Monitoring Modern Power Grids with  Real-Time Dynamic Harmonic Estimation\n\nInput:\n\nsol::AbstractDTFTSolution   | DTFT solution struct [-]\nH::Int                      | harmonic number [-], default=1\n\nOutput:\n\nf::Vector{<:Real}           | frequency fₕ(t) [Hz]\n\n\n\n\n\n","category":"function"},{"location":"util/#Rate-Of-Change-Of-Frequency-(ROCOF)","page":"Utilities","title":"Rate-Of-Change-Of-Frequency (ROCOF)","text":"","category":"section"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.r(sol::TFT.AbstractDTFTSolution, H::Int=1)","category":"page"},{"location":"util/#TFT.r","page":"Utilities","title":"TFT.r","text":"TFT.r(sol::TFT.AbstractDTFTSolution, H::Int=1)\n\nShorthand function to obtain the rate-of-change-of-frequency of the  zeroth-degree derivative of the Hth-harmonic phasor, dispatching to  rocof(sol,H).\n\n\n\n\n\n","category":"function"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.rocof(sol::TFT.AbstractDTFTSolution, H::Int=1)","category":"page"},{"location":"util/#TFT.rocof","page":"Utilities","title":"TFT.rocof","text":"TFT.rocof(sol::TFT.AbstractDTFTSolution, H::Int=1)\n\nFunction to obtain the rate-of-change-of-frequency of the zeroth-degree  derivative of the Hth-harmonic phasor.\n\nrₕ(t) = ϕₕ⁽²⁾(t) / (2 π)² ∈ 𝐑\n\nSee: Fast Taylor-Fourier Transform for Monitoring Modern Power Grids with  Real-Time Dynamic Harmonic Estimation\n\nInput:\n\nsol::AbstractDTFTSolution   | DTFT solution struct [-]\nH::Int                      | harmonic number [-], default=1\n\nOutput:\n\nr::Vector{<:Real}           | rocof rₕ(t) [Hz/s]\n\n\n\n\n\n","category":"function"},{"location":"util/#Dynamic-Phasor","page":"Utilities","title":"Dynamic Phasor","text":"","category":"section"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.ξ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.ξ","page":"Utilities","title":"TFT.ξ","text":"TFT.ξ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nShorthand function to obtain the Dth-degree derivative of the Hth-harmonic  dynamic phasor, dispatching to phasor(sol,D,H).\n\n\n\n\n\n","category":"function"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.phasor(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.phasor","page":"Utilities","title":"TFT.phasor","text":"TFT.phasor(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nFunction to obtain the D-th-degree derivative of the Hth-harmonic dynamic  phasor. For the zeroth-harmonic dynamic phasor, only the real part of the  dynamic phasor is returned.\n\n∀ h ∈ {0}:\n    ξ₀⁽ᴰ⁾(t) ∈ ℝ\n∀ h ∈ 𝓗/{0}\n    ξₕ⁽ᴰ⁾(t) ∈ ℂ\n\nInput:\n\nsol::AbstractDTFTSolution   | DTFT solution struct [-]\nD::Int                      | degree of the derivative [-], default=0\nH::Int                      | harmonic number [-], default=1\n\nOutput:\n\nξ::Vector{<:Complex}        | dynamic phasor ξₕ⁽ᴰ⁾(t) [?]\n\n\n\n\n\n","category":"function"},{"location":"util/#Anti-Rotating-Dynamic-Phasor","page":"Utilities","title":"Anti-Rotating Dynamic Phasor","text":"","category":"section"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.ψ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.ψ","page":"Utilities","title":"TFT.ψ","text":"TFT.ψ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nShorthand function to obtain the Dth-degree derivative of the Hth-harmonic anti- rotating dynamic phasor, dispatching to ar_phasor(sol,D,H).\n\n\n\n\n\n","category":"function"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.ar_phasor(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.ar_phasor","page":"Utilities","title":"TFT.ar_phasor","text":"TFT.ar_phasor(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nFunction to obtain the D-th-degree derivative of the Hth-harmonic anti-rotating  dynamic phasor.\n\n∀ h ∈ {0}:\n    ψ₀⁽ᴰ⁾(t) = ξ₀⁽ᴰ⁾(t) ∈ ℝ\n∀ h ∈ 𝓗/{0}:\n    ψₕ⁽ᴰ⁾(t) = ξₕ⁽ᴰ⁾(t) exp(-im ωₕ t) ∈ 𝐂\n\nInput:\n\nsol::AbstractDTFTSolution   | DTFT solution struct [-]\nD::Int                      | degree of the derivative [-], default=0\nH::Int                      | harmonic number [-], default=1\n\nOutput:\n\nψ::Vector{<:Complex}        | anti-rotating dynamic phasor ψₕ⁽ᴰ⁾(t) [?]\n\n\n\n\n\n","category":"function"},{"location":"util/#Signal","page":"Utilities","title":"Signal","text":"","category":"section"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.signal(sol::TFT.AbstractDTFTSolution)","category":"page"},{"location":"util/#TFT.signal-Tuple{TFT.AbstractDTFTSolution}","page":"Utilities","title":"TFT.signal","text":"TFT.signal(sol::TFT.AbstractDTFTSolution)\n\nFunction to obtain the overall signal.\n\ns(t) = ∑ₕ ξₕ⁽⁰⁾(t) + conj(ξₕ⁽⁰⁾(t)) ∈ 𝐑\n\nInput:\n\nsol::AbstractDTFTSolution   | DTFT solution struct [-]\n\nOutput:\n\ns::Vector{<:Real}           | signal s(t) [?]\n\n\n\n\n\n","category":"method"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.signal(sol::TFT.AbstractDTFTSolution, H::Int)","category":"page"},{"location":"util/#TFT.signal-Tuple{TFT.AbstractDTFTSolution, Int64}","page":"Utilities","title":"TFT.signal","text":"TFT.signal(sol::TFT.AbstractDTFTSolution, H::Int)\n\nFunction to obtain the Hth-harmonic signal. \n\n∀ h ∈ {0}:\n    s₀(t) = ξ₀⁽⁰⁾(t) ∈ 𝐑\n∀ h ∈ H/{0}: \n    sₕ(t) = ℜ[ξₕ⁽⁰⁾(t) + conj(ξₕ⁽⁰⁾(t))] ∈ 𝐑\n\nNote: In its essence, the real operator ℜ[]is unnecessary, however, it is  used to convert the complex numberx + j0to a real numberx`. \n\nInput:\n\nsol::AbstractDTFTSolution   | DTFT solution struct [-]\nH::Int                      | harmonic number [-]\n\nOutput:\n\ns::Vector{<:Real}           | signal sₕ(t) [?]\n\n\n\n\n\n","category":"method"},{"location":"util/#Error","page":"Utilities","title":"Error","text":"","category":"section"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.error(sol::TFT.AbstractDTFTSolution)","category":"page"},{"location":"util/#TFT.error-Tuple{TFT.AbstractDTFTSolution}","page":"Utilities","title":"TFT.error","text":"TFT.error(sol::TFT.AbstractDTFTSolution)\n\nFunction to obtain the error between the input and computed signal.\n\nInput:\n\nsol::AbstractDTFTSolution   | DTFT solution struct [-]\n\nOutput:\n\ne::Vector{<:Real}           | error [?]\n\n\n\n\n\n","category":"method"}]
}
