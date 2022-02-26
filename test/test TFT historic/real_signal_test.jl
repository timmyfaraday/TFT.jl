# using pkgs
using CSV 
using DataFrames
using Plots
using TFT

# input parameters
N   = 1000
F   = 50.0
K   = 9
H   = collect(1:20)

# read-in data
df = DataFrame(CSV.File("C:/Users/VANACKT5/OneDrive - BASF/Work/Power Harmonics/Data/INRUSH_NT2_D0542.CSV", header=14))

# clean DataFrame
rename!(df,[:Time, :Ia, :Ib, :Ic, :Van, :Vbn, :Vcn, :Vs])

# time idx 
idx_fr = findlast(x -> x .<= 0.00, df.Time)
idx_to = findlast(x -> x .<= 1.20, df.Time)
idx = idx_fr:idx_to
t = df.Time[idx]

# initialize a dictionary for the signals
sym = [:Ia, :Ib, :Ic, :Van, :Vbn, :Vcn]
signal = Dict(ns => Dict{Int,Vector{<:Complex}}() for ns in sym)

# 1) determine all amplitudes and angles of Ia, Ib, Ic, Van, Vbn, Vcn 
# 2) reconstruct time signal: a⁰ .* exp(im .* φ⁰)
for ns in sym, nh in H
    a⁰, ϕ⁰, a¹, ϕ¹, a², ϕ² = harmonic_estimator(df[!,ns], K, N, nh, F)
    signal[ns][nh] = a⁰[idx] .* exp.(im .* ϕ⁰[idx])
end

# 3) determine the sequence components through fortescue
a = exp(2.0 * pi * im / 3.0)
for ns in [:V0, :V1, :V2, :I0, :I1, :I2]
    signal[ns] = Dict{Int,Vector{<:Complex}}()
end
for nh in H
    signal[:V0][nh] = signal[:Van][nh] .+ signal[:Vbn][nh] .+ signal[:Vcn][nh]
    signal[:V1][nh] = signal[:Van][nh] .+ a .* signal[:Vbn][nh] .+ a^2 .* signal[:Vcn][nh]
    signal[:V2][nh] = signal[:Van][nh] .+ a^2 .* signal[:Vbn][nh] .+ a .* signal[:Vcn][nh]
end
for nh in H
    signal[:I0][nh] = signal[:Ia][nh] .+ signal[:Ib][nh] .+ signal[:Ic][nh]
    signal[:I1][nh] = signal[:Ia][nh] .+ a .* signal[:Ib][nh] .+ a^2 .* signal[:Ic][nh]
    signal[:I2][nh] = signal[:Ia][nh] .+ a^2 .* signal[:Ib][nh] .+ a .* signal[:Ic][nh]
end

T = t[findlast(x -> x .<= 0.25, t):findlast(x -> x .<= 1.00, t)]

nh = 10
nd = 9

pol_V = fit(t, signal[:V1][nh], nd)
pol_I = fit(t, signal[:I1][nh], nd)

dp_V = derivative(pol_V)
dp_I = derivative(pol_I)

plot(t, abs.(dp_V.(t) ./ dp_I.(t)),
        #label=["base" "derivative"],
        title="Positive sequence $(nh)th-harmonic network impedance \n at primary of D0540-NT2 (polynomial degree = $(nd))",
        xlabel="Time [s]",
        ylabel="|Z₁| [Ω]",
        ylim=(0.0,50.0))


# Procedure
# 1) determine all amplitudes and angles of Ia, Ib, Ic, Van, Vbn, Vcn 
# 2) reconstruct time signal: a⁰ .* exp(im .* φ⁰)
# 3) determine the sequence components through fortescue
#       X0 = 1/3 .* (Xa .+ Xb .+ Xc)
#       X1 = 1/3 .* (Xa .+ a.*Xb .+ a^2.*Xc)
#       X2 = 1/3 .* (Xa .+ a^2.*Xb .+ a.*Xc)
# 4) fit a polynomial through all sequence components Xs
# 5) determine the derivative of all sequence components dXs
# 6) determine the sequence impedance Zs = dVs / dIs 

# notes to RL 
# 1) one second before switch-in, nine seconds after switch-in
# 2) play around with the sampling rate
# 3) explenation of the multipliers