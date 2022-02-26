using TFT
using Plots

# input parameters
N   = 1000
F   = 50.0
K   = 9
H   = 1

# additional parameters
ω   = 2 * pi * F

# time range
T   = 1 / F
ΔT  = T / N
t   = 0.0:ΔT:1.5

# analytical signal and derivative
s(t)  = (10 - t^2 + t^3) * cos(ω * t + pi/2 * t^2)
ds(t) = (3 * t^2 - 2 * t) * cos(ω * t + pi/2 * t^2) - (ω + pi *t) * (10 - t^2 + t^3) * sin(ω * t + pi/2 * t^2)

# TFT
ξ⁰, ξ¹, ξ², a⁰, ϕ⁰, a¹, ϕ¹, a², ϕ² = harmonic_estimator(s.(t), K, N, H, F)

# discrete signal and derivative reconstruction
S   = real(a⁰ .* exp.(im .* ϕ⁰) .* exp.(im .* ω .* t))
dS  = real((a¹ .+ a⁰ .* im .* ϕ¹ .+ a⁰ .* im .* ω) .* exp.(im * ϕ⁰) .* exp.(im .* ω .* t))

# discrete phasor
P   = a⁰ .* exp.(im .* ϕ⁰)
dP  = a¹ .* exp.(im .* ϕ⁰) .+ a⁰ .* im .* ϕ¹ .* exp.(im .* ϕ⁰)

# signal plots
plot(t, [s.(t), S, a⁰])
plot(t, [ds.(t), dS, a¹ .+ a⁰ .* ϕ¹ .+ a⁰ .* ω])

# errors
Δ   = s.(t) - S
dΔ  = ds.(t) - dS

# amplitude and phase plots
τ   = 0.1:ΔT:1.4
idx = 5000:70000

# error plots 
plot(τ, Δ[idx])
plot(τ, dΔ[idx])

plot(τ, [a⁰[idx], 10 .- τ.^2 + τ.^3])
plot(τ, [ϕ⁰[idx], pi / 2 .* τ.^2])

plot(τ, [a¹[idx], -2.0 .* τ .+ 3.0 .* τ.^2])
plot(τ, [ϕ¹[idx], pi .* τ])

plot(τ, [a²[idx], -2.0 .+ 6.0 .* τ])
plot(τ, [ϕ²[idx], pi .* ones(length(τ))])