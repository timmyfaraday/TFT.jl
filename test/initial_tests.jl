# using pkgs
using TFT
using Plots

# input parameters
N   = 100
F   = 50.0
K   = 9
H   = 7

# build sampling time
T   = 1 / F
Δt  = T / N 
t   = 0.0:Δt:3.0-Δt

# build the signal `s`
s₁  = 10.0 .* sin.((2 * pi * F) .* t)
s₇  = (10.0 .- ifelse.(0.2 .<= t .<= 0.4,3.0,0.0)) .* sin.((2 * pi * F * 7) .* t)
s  = s₁ .+ s₇

# get the info of the harmonic `H`
a⁰, ϕ⁰, a¹, ϕ¹, a², ϕ² = harmonic_state_estimator(s, K, N, H, F)
