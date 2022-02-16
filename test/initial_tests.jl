# using pkgs
using TFT
using Plots

# input parameters
N   = 100
F   = 50.0
K   = 9
H   = 1

# build sampling time
T   = 1 / F
Δt  = T / N 
t   = 0.0:Δt:1.0-Δt

# build the signal `s`
s   = 10.0 .* sin.((2 * pi * F) .* t)

# get the info of the harmonic `H`
a⁰, ϕ⁰, a¹, ϕ¹, a², ϕ² = harmonic_state_estimator(s, K, N, H, F)

# plot 
plot(t, s)
plot!(t, a⁰)