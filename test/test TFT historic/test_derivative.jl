t = 0.0:0.001:1.0
ω = 2 * π * 50.0

O(t) = (10 - t) * cos(ω * t + pi/2 * t)
S(t) = (10 - t) * exp(im * pi/2 * t) * exp(im * ω * t)

dO(t) = -(ω + pi/2) * (10 - t) * sin(ω * t + pi * t / 2) - cos(ω * t + pi * t / 2)
dS(t) = (-exp(im * pi/2 * t) + im * pi/2 * (10 - t) * exp(im * pi/2 * t) + im * ω * (10 - t) * exp(im * pi/2 * t)) * exp(im * ω * t)

using TFT
using Plots

s = O.(t)

ξ⁰, ξ¹, ξ², a⁰, ϕ⁰, a¹, ϕ¹, a², ϕ² = harmonic_estimator(s, 9, 20, 1, 50.0)

dV = (ξ¹ .+ im .* ω .* ξ⁰) .* exp.(im .* ω .* t)

plot(t,real.(dV))
plot!(t,real(dS.(t)))


v(t) = cos(ω * t)
dv(t) = - ω * sin(ω * t)

ξ⁰, ξ¹, ξ², a⁰, ϕ⁰, a¹, ϕ¹, a², ϕ² = harmonic_estimator(v.(t), 9, 20, 1, 50.0)

V = real(a⁰ .* exp.(im * ϕ⁰) .* exp.(im .* ω .* t))
dV = real((a¹ .* exp.(im * ϕ⁰) .+ a⁰ .* ϕ¹ .* exp.(im * ϕ⁰) .+ im .* ω .* a⁰ .* exp.(im * ϕ⁰)) .* exp.(im .* ω .* t))

### this is the equation!!!