# using pkgs 
using CSV
using DataFrames
using Plots 
using TFT

# read-in data
df = DataFrame(CSV.File("C:/Users/VANACKT5/OneDrive - BASF/Work/Power Harmonics/Data/INRUSH_NT2_D0542.CSV", header=14))

# clean DataFrame
rename!(df,[:Time, :Ia, :Ib, :Ic, :Van, :Vbn, :Vcn, :Vs])

# relevant signals
sym = [:Ia, :Ib, :Ic, :Van, :Vbn, :Vcn]

# perform the TFT
sol = Dict(ns => tft(df[!,ns], df.Time, collect(1:20), 1, 50.0, 9) for ns in sym)

# find the phasors of the first derivative
pha = Dict(ns => Dict(nh => phasor(sol[ns], nh, 1) for nh in sol[ns].prob.h) for ns in sym)

# find the sequence components of the phasors
a = exp.(im .* 2.0 .* pi ./ 3.0)
seq = Dict( :I0 => Dict(nh => pha[:Ia][nh] .+ pha[:Ib][nh] + pha[:Ic][nh] for nh in sol.prob.h),
            :I1 => Dict(nh => pha[:Ia][nh] .+ a^2 .* pha[:Ib][nh] + a .* pha[:Ic][nh] for nh in sol.prob.h),
            :I2 => Dict(nh => pha[:Ia][nh] .+ a .* pha[:Ib][nh] + a^2 .* pha[:Ic][nh] for nh in sol.prob.h),
            :V0 => Dict(nh => pha[:Van][nh] .+ pha[:Vbn][nh] + pha[:Vcn][nh] for nh in sol.prob.h),
            :V1 => Dict(nh => pha[:Van][nh] .+ a^2 .* pha[:Vbn][nh] + a .* pha[:Vcn][nh] for nh in sol.prob.h),
            :V2 => Dict(nh => pha[:Van][nh] .+ a .* pha[:Vbn][nh] + a^2 .* pha[:Vcn][nh] for nh in sol.prob.h))

# determine the sequence impedances
imp = Dict( :Z0 => Dict(nh => seq[:V0][nh] ./ seq[:I0][nh] for nh in sol.prob.h),
            :Z1 => Dict(nh => seq[:V1][nh] ./ seq[:I1][nh] for nh in sol.prob.h),
            :Z2 => Dict(nh => seq[:V2][nh] ./ seq[:I2][nh] for nh in sol.prob.h))

# determine the phase impedances
fin = Dict( :Za => Dict(nh => imp[:Z0][nh] .+ imp[:Z1][nh] + imp[:Z2][nh] for nh in sol.prob.h),
            :Zb => Dict(nh => imp[:Z0][nh] .+ a .* imp[:Z1][nh] + a^2 .* imp[:Z2][nh] for nh in sol.prob.h),
            :Zc => Dict(nh => imp[:Z0][nh] .+ a^2 .* imp[:Z1][nh] + a .* imp[:Z2][nh] for nh in sol.prob.h))
