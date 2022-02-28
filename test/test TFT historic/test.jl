using TFT 
using Plots

# 
A(t)    = 10.0 - t
dA(t)   = -t

Φ(t)    = pi/2
dΦ(t)   = 0.0

F   = 50.0
ω   = 2.0 * pi * F 

D   = 1
K   = 9

t   = -0.1:0.0001:1.1

S(t)   = real(A.(t) .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t))
dS(t)  = real(  dA.(t) .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t) .+
                A.(t) .* im .* dΦ.(t) .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t) .+
                A.(t) .* exp.(im .* Φ.(t)) .* im .* ω .* exp.(im .* ω .* t))

sol = tft(S.(t), collect(t), [1], 1, F, K)

A0 = amplitude(sol,0,1)
A1 = amplitude(sol,1,1)

Φ0 = TFT.angle(sol,0,1)
Φ1 = TFT.angle(sol,1,1)

P0 = phasor(sol,0,1)
P1 = phasor(sol,1,1)

S0 = signal(sol,0,1)
S1 = signal(sol,1,1)