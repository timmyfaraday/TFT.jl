################################################################################
#  Copyright 2022, Tom Van Acker                                               #
################################################################################
# TFT.jl                                                                       #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TFT.jl                                    #
################################################################################

"""
    TFT.signal(TFT.AbstractDTFTSolution, deg::Int, har)

xxx
"""
function signal(sol::AbstractDTFTSolution, D::Int=0, H::Int=1)
    D == 0 && return phasor(sol,D,H) .* exp.(im .* 2.0 .* pi .* H .* sol.prob.F .* sol.t)
end

"""
    TFT.phasor(TFT.AbstractDFTFSolution, deg::Int, har::Int)

xxx
"""
function phasor(sol::AbstractDFTFSolution, D::Int=0, H::Int=1)
    D == 0 && return amplitude(sol,D,H) .* exp.(im .* angle(sol,D,H))
end

"""
    TFT.amplitude(TFT.AbstractDFTFSolution, deg::Int, har::Int)

xxx
"""
function amplitude(sol::AbstractDFTFSolution, D::Int=0, H::Int=1)
    D == 0 && return 2.0 .* abs.(sol.x[0,H,:])
    D == 1 && return 2.0 .* real.(sol.x[1,H,:])
    D == 2 && return 2.0 .* (real.(sol.x[1,H,:] .+ amplitude(sol,1,H) .* angle(sol,1,H)))
end

"""
    TFT.angle(TFT.AbstractDFTFSolution, deg::Int, har::Int)

xxx
"""
function angle(sol::AbstractDFTFSolution, D::Int=0, H::Int=1)
    D == 0 && return angle.()
    D == 1 && return imag.(sol.x[1,H,:]) ./ amplitude(sol,0,H)
    D == 2 && return (imag.(sol.x[2,H,:]) .- 2.0 .* amplitude(sol,1,H) .* angle(sol,1,H)) ./ amplitude(sol,0,H)
end

"""
    TFT.frequency(TFT.AbstractDTFTSolution, har::Int)

xxx
"""
function frequency(sol::AbstractDFTFSolution, H::Int=1)
    return H * sol.prob.F .+ angle(sol,1,H) ./ (2.0 * pi)
end

"""
    TFT.rocof(sol::TFT.AbstractDFTFSolution, har::Int)

xxx
"""
function rocof(sol::AbstractDFTFSolution, H::Int=1)
    return angle(sol,2,H) ./ (2.0 * pi)
end


"""
    TFT.error(sol::TFT.AbstractDFTFSolution)

xxx
"""
function error(sol::AbstractDFTFSolution)

end