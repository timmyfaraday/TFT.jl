var documenterSearchIndex = {"docs":
[{"location":"math/#Mathematical-Background-on-Taylor-Fourier-Transform","page":"Mathematical Background","title":"Mathematical Background on Taylor-Fourier Transform","text":"","category":"section"},{"location":"math/#Nomenclature","page":"Mathematical Background","title":"Nomenclature","text":"","category":"section"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"Symbol Description\ns signal\na amplitude\nϕ phase [(rad)]\nφ anti-rotating phase [(rad)]\nξ dynamic phasor\nψ anti-rotating dynamic phasor","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"Set Description\nh ∈ H set of harmonic number [Int]","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"Parameter Description\nF frequency [Hz]\nω angular frequency [(rad)/s]","category":"page"},{"location":"math/#Taylor-Fourier-Transform","page":"Mathematical Background","title":"Taylor-Fourier Transform","text":"","category":"section"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"The Taylor-Fourier Transform function tft() gives the up-to-Dth derivative of  the hth-harmonic dynamic phasors ξₕ⁽ᵈ⁾(t), ∀ d ∈ {0,..,D}, h ∈ H. ","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"The zeroth derivative of the hth-harmonic dynamic phasor ξₕ⁽⁰⁾(t) is given by,","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"beginaligned\n    xi^(0)_h(t)        = fraca^(0)_h(t)2  exp(j phi^(0)_h(t)) \n                            = fraca^(0)_h(t)2  exp(j varphi^(0)_h(t))  exp(j omega_h t) \n                            = psi^(0)_h(t)  exp(j omega_h t)\n\nendaligned","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"where ψₕ⁽⁰⁾(t) denotes the anti-rotating zeroth derivative of the hth-harmonic  dynamic phasor. ","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"The corresponding amplitude aₕ⁽⁰⁾(t), phase ϕₕ⁽⁰⁾(t) and anti-rotating phase  φₕ⁽⁰⁾(t) are, respectively,","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"beginaligned\n    a^(0)_h(t)          = 2 cdot xi^(0)_h(t) \n    phi^(0)_h(t)       = angle big( xi^(0)_h(t) big) \n    varphi^(0)_h(t)    = angle big( xi^(0)_h(t)  exp(-j omega_h t) big) = angle big( psi^(0)_h(t) big)\nendaligned","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"The corresponding hth-harmonic signal sₕ⁽⁰⁾(t) and overall signal s⁽⁰⁾(t) are,  respectively,","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"beginaligned\n    s^(0)_h(t)          = xi^(0)_h(t) + textconjbig( xi^(0)_h(t) big) \n    s^(0)(t)              = sum_h in H xi^(0)_h(t) + textconjbig( xi^(0)_h(t) big)\nendaligned","category":"page"},{"location":"math/","page":"Mathematical Background","title":"Mathematical Background","text":"where conj(ξₕ⁽⁰⁾(t)) denotes the complex conjugate of the zeroth derivative of the  hth-harmonic dynamic phasor.","category":"page"},{"location":"#TFT.jl","page":"Home","title":"TFT.jl","text":"","category":"section"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"TFT.jl is a Julia package for Taylor-Fourier Transform.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The latest stable release of TFT.jl can be installed using the Julia package  manager:","category":"page"},{"location":"","page":"Home","title":"Home","text":"] add \"https://github.com/timmyfaraday/TFT.jl\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"In order to test whether the package works, run:","category":"page"},{"location":"","page":"Home","title":"Home","text":"] test TFT","category":"page"},{"location":"util/#Utilities","page":"Utilities","title":"Utilities","text":"","category":"section"},{"location":"util/","page":"Utilities","title":"Utilities","text":"A number of functions are made available to the user to retrieve specific components of the TFT solution","category":"page"},{"location":"util/#Quick-Links","page":"Utilities","title":"Quick Links","text":"","category":"section"},{"location":"util/","page":"Utilities","title":"Utilities","text":"Utilities\nQuick Links\nAmplitude\nPhase\nAnti-Rotating Phase\nDynamic Phasor\nAnti-Rotating Dynamic Phasor\nSignal\nError","category":"page"},{"location":"util/#Amplitude","page":"Utilities","title":"Amplitude","text":"","category":"section"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.a(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.a","page":"Utilities","title":"TFT.a","text":"TFT.a(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nShorthand function to obtain the amplitude of the Dth-degree derivative of the Hth-harmonic phasor, dispatching to amplitude(sol,D,H).\n\n\n\n\n\n","category":"function"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.amplitude(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.amplitude","page":"Utilities","title":"TFT.amplitude","text":"TFT.amplitude(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nFunction to obtain the amplitude of the Dth-degree derivative of the  Hth-harmonic phasor.\n\naₕ⁽⁰⁾(t) = 2 ⋅ |ξₕ⁽⁰⁾(t)|\n# aₕ⁽¹⁾(t) = 2 ⋅ ℜ[ξₕ⁽¹⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] \n# aₕ⁽²⁾(t) = 2 ⋅ ℜ[ξₕ⁽²⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] + aₕ⁽⁰⁾(t) ⋅ φₕ⁽¹⁾(t)^2\n\nSee: Assessing Synchrophasor Estimates of an Event Captured by a Phasor  Measurement Unit, pg. 3112\n\nInput:\n\nsol::AbstractDTFTSolution   | DTFT solution struct [-]\nD::Int                      | degree of the derivative [-], default=0\nH::Int                      | harmonic number [-], default=1\n\nOutput:\n\na::Vector{<:Real}           | amplitude aₕ⁽ᴰ⁾(t) [?]\n\n\n\n\n\n","category":"function"},{"location":"util/#Phase","page":"Utilities","title":"Phase","text":"","category":"section"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.ϕ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.ϕ","page":"Utilities","title":"TFT.ϕ","text":"TFT.ϕ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nShorthand function to obtain the phase of the Dth-degree derivative  of the Hth-harmonic phasor, dispatching to phase(sol, D, H).\n\n\n\n\n\n","category":"function"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.phase(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.phase","page":"Utilities","title":"TFT.phase","text":"TFT.phase(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nFunction to obtain the alternative angle of the Dth-degree derivative of the  Hth-harmonic phasor.\n\nϕₕ⁽⁰⁾(t) = ∠[pₕ⁽⁰⁾(t)]\n\nSee: Assessing Synchrophasor Estimates of an Event Captured by a Phasor  Measurement Unit, pg. 3112\n\nInput:\n\nsol::AbstractDTFTSolution   | DTFT solution struct [-]\nH::Int                      | harmonic number [-], default=1\n\nOutput:\n\nϕ::Vector{<:Real}           | phase ϕₕ⁽ᴰ⁾(t) [(rad)]\n\n\n\n\n\n","category":"function"},{"location":"util/#Anti-Rotating-Phase","page":"Utilities","title":"Anti-Rotating Phase","text":"","category":"section"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.φ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.φ","page":"Utilities","title":"TFT.φ","text":"TFT.φ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nShorthand function to obtain the anti-rotating phase of the Dth-degree  derivative of the Hth-harmonic phasor, dispatching to angle(sol,D,H).\n\n\n\n\n\n","category":"function"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.ar_phase(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.ar_phase","page":"Utilities","title":"TFT.ar_phase","text":"TFT.ar_phase(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nFunction to obtain the anti-rotating phase of the Dth-degree derivative of the  Hth-harmonic phasor.\n\nφₕ⁽⁰⁾(t) = ∠[ξₕ⁽⁰⁾(t) ⋅ exp(-2.0 ⋅ im ⋅ pi ⋅ H / N) ⋅ 0:length(s(t))-1]\n# φₕ⁽¹⁾(t) = ℑ[pₕ⁽¹⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] / aₕ⁽⁰⁾(t)\n# φₕ⁽²⁾(t) = (ℑ[pₕ⁽²⁾(t) ⋅ exp(-im ⋅ ϕₕ⁽⁰⁾(t))] + 2 ⋅ aₕ⁽¹⁾(t) ⋅ φₕ⁽¹⁾(t)) / aₕ⁽⁰⁾(t)\n\nSee: Assessing Synchrophasor Estimates of an Event Captured by a Phasor  Measurement Unit, pg. 3112\n\nInput:\n\nsol::AbstractDTFTSolution   | DTFT solution struct [-]\nD::Int                      | degree of the derivative [-], default=0\nH::Int                      | harmonic number [-], default=1\n\nOutput:\n\nφ::Vector{<:Real}           | anti-rotating phase φₕ⁽ᴰ⁾(t) [(rad)]\n\n\n\n\n\n","category":"function"},{"location":"util/#Dynamic-Phasor","page":"Utilities","title":"Dynamic Phasor","text":"","category":"section"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.ξ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.ξ","page":"Utilities","title":"TFT.ξ","text":"TFT.ξ(TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nShorthand function to obtain the Dth-degree derivative of the Hth-harmonic  dynamic phasor, dispatching to phasor(sol,D,H).\n\n\n\n\n\n","category":"function"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.phasor(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.phasor","page":"Utilities","title":"TFT.phasor","text":"TFT.phasor(TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nFunction to obtain the D-th-degree derivative of the Hth-harmonic dynamic  phasor.\n\nInput:\n\nsol::AbstractDTFTSolution   | DTFT solution struct [-]\nD::Int                      | degree of the derivative [-], default=0\nH::Int                      | harmonic number [-], default=1\n\nOutput:\n\nphasor::Vector{<:Complex}     | dynamic phasor ξₕ⁽ᴰ⁾(t) [?]\n\n\n\n\n\n","category":"function"},{"location":"util/#Anti-Rotating-Dynamic-Phasor","page":"Utilities","title":"Anti-Rotating Dynamic Phasor","text":"","category":"section"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.ψ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.ψ","page":"Utilities","title":"TFT.ψ","text":"TFT.ψ(TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nShorthand function to obtain the Dth-degree derivative of the Hth-harmonic anti- rotating dynamic phasor, dispatching to ar_phasor(sol,D,H).\n\n\n\n\n\n","category":"function"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.ar_phasor(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)","category":"page"},{"location":"util/#TFT.ar_phasor","page":"Utilities","title":"TFT.ar_phasor","text":"TFT.ar_phasor(TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)\n\nFunction to obtain the D-th-degree derivative of the Hth-harmonic anti-rotating  dynamic phasor.\n\nψₕ⁽ᴰ⁾(t) = ξₕ⁽ᴰ⁾(t) exp(-im ωₕ t)\n\nInput:\n\nsol::AbstractDTFTSolution   | DTFT solution struct [-]\nD::Int                      | degree of the derivative [-], default=0\nH::Int                      | harmonic number [-], default=1\n\nOutput:\n\nar_phasor::Vector{<:Complex}  | anti-rotating dynamic phasor ψₕ⁽ᴰ⁾(t) [?]\n\n\n\n\n\n","category":"function"},{"location":"util/#Signal","page":"Utilities","title":"Signal","text":"","category":"section"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.signal(sol::TFT.AbstractDTFTSolution)","category":"page"},{"location":"util/#TFT.signal-Tuple{TFT.AbstractDTFTSolution}","page":"Utilities","title":"TFT.signal","text":"TFT.signal(TFT.AbstractDTFTSolution)\n\nFunction to obtain the overall signal.\n\ns(t) = ∑ₕ ξₕ⁽⁰⁾(t) + conj(ξₕ⁽⁰⁾(t))\n\nInput:\n\nsol::AbstractDTFTSolution   | DTFT solution struct [-]\n\nOutput:\n\nsignal::Vector{<:Real}      | signal s(t) [?]\n\n\n\n\n\n","category":"method"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.signal(sol::TFT.AbstractDTFTSolution, H::Int)","category":"page"},{"location":"util/#TFT.signal-Tuple{TFT.AbstractDTFTSolution, Int64}","page":"Utilities","title":"TFT.signal","text":"TFT.signal(TFT.AbstractDTFTSolution, H::Int)\n\nFunction to obtain the Hth-harmonic signal. \n\nsₕ(t) = ξₕ⁽⁰⁾(t) + conj(ξₕ⁽⁰⁾(t))\n\nInput:\n\nsol::AbstractDTFTSolution   | DTFT solution struct [-]\nH::Int                      | harmonic number [-]\n\nOutput:\n\nsignal::Vector{<:Real}      | signal sₕ(t) [?]\n\n\n\n\n\n","category":"method"},{"location":"util/#Error","page":"Utilities","title":"Error","text":"","category":"section"},{"location":"util/","page":"Utilities","title":"Utilities","text":"TFT.error(sol::TFT.AbstractDTFTSolution)","category":"page"},{"location":"util/#TFT.error-Tuple{TFT.AbstractDTFTSolution}","page":"Utilities","title":"TFT.error","text":"TFT.error(sol::TFT.AbstractDTFTSolution)\n\nFunction to obtain the error between the input and computed signal.\n\nInput:\n\nsol::AbstractDTFTSolution   | DTFT solution struct [-]\n\nOutput:\n\ne::Vector{<:Real}           | error [?]\n\n\n\n\n\n","category":"method"}]
}
