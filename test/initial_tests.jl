# using pkgs
using TFT
using Plots

# input parameters
N   = 100
F   = 50.0
K   = 9
H   = [1,7]

# build sampling time
T   = 1 / F
Δt  = T / N 
t   = 0.0:Δt:1.0-Δt

# build the signal `s`
s₁  = 10.0 .* sin.((2 * pi * F) .* t)
s₇  = (3.0 .- ifelse.(0.2 .<= t .<= 0.4,2.0,0.0)) .* sin.((2 * pi * F * 7) .* t .+ t)
s  = s₁ .+ s₇

# get the info of the harmonic `H`
a⁰, ϕ⁰  = Dict{Int,Vector{<:Real}}(), Dict{Int,Vector{<:Real}}()
a¹, ϕ¹ = Dict{Int,Vector{<:Real}}(), Dict{Int,Vector{<:Real}}()
a², ϕ² = Dict{Int,Vector{<:Real}}(), Dict{Int,Vector{<:Real}}()
@time for h in H
    a⁰[h], ϕ⁰[h], a¹[h], ϕ¹[h], a²[h], ϕ²[h] = harmonic_estimator(s, K, N, h, F)
end


plot(t,s₁)
plot!(t,a⁰[1])