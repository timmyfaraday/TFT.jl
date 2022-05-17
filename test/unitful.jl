################################################################################
#  Copyright 2022, Tom Van Acker (BASF Antwerp)                                #
################################################################################
# TFT.jl                                                                       #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TFT.jl                                    #
################################################################################

@testset "Unitful" begin
    
    # fundamental frequency and angular frequency
    F   = 50.0u"Hz"
    ω   = 2 * pi * F

    # tft input
    D   = 2
    K   = 9

    # discrete time
    t   = (0.0:0.001:1.0)u"s"
    tm  = 0.536u"s"
    idm = findfirst(x -> x == tm, t)

    # tolerances
    atol_f  = (atol)u"Hz"
    atol_r  = (atol)u"Hz^2"

    @testset "Fundamental Aperiodic Linear Signal" begin 
        # tolerances
        atol_u  = (atol)u"V"

        # input - amplitude and anti-rotating phase
        A(t)    = 10.0u"V" .- (t)u"V/s"
        Φ(t)    = pi / 2 .* (t)u"1/s"

        # derivatives - amplitude and anti-rotating phase
        dA(t)   = -1.0u"V"
        d2A(t)  = 0.0u"V"
        dΦ(t)   = pi / 2
        d2Φ(t)  = 0.0

        # derived input
        Φr(t)   = rem.(ω .*t .+ Φ.(t) .+ pi, 2 * pi) .- pi
        Fr(t)   = F + dΦ(t)u"Hz" / (2 * pi)
        Rr(t)   = d2Φ(t)u"Hz^2" / (2 * pi)^2
        Ξ(t)    = A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t)
        Ψ(t)    = A.(t) ./ 2 .* exp.(im .* Φ.(t))
        S(t)    = real.(A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t) .+ 
                  conj.(A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t)))

        # perform taylor-fourier transform
        sol = tft(S.(t), collect(t), [1], D, F, K)

        # tests
        ## amplitude
        @test isapprox(TFT.a(sol,0,1)[idm], A(tm), atol=atol_u)
        @test isapprox(TFT.a(sol,1,1)[idm], dA(tm), atol=atol_u)
        @test isapprox(TFT.a(sol,2,1)[idm], d2A(tm), atol=atol_u)
        ## phase
        @test isapprox(TFT.ϕ(sol,0,1)[idm], Φr(tm), atol=atol)
        @test isapprox(TFT.ϕ(sol,1,1)[idm], dΦ(tm), atol=atol)
        @test isapprox(TFT.ϕ(sol,2,1)[idm], d2Φ(tm), atol=atol)
        ## anti-rotating phase
        @test isapprox(TFT.φ(sol,0,1)[idm], Φ(tm), atol=atol)
        @test isapprox(TFT.φ(sol,1,1)[idm], dΦ(tm), atol=atol)
        @test isapprox(TFT.φ(sol,2,1)[idm], d2Φ(tm), atol=atol)
        ## frequency 
        @test isapprox(TFT.f(sol,1)[idm], Fr(tm), atol=atol_f)
        ## rocof
        @test isapprox(TFT.r(sol,1)[idm], Rr(tm), atol=atol_r)
        ## dynamic phasor
        @test isapprox(TFT.ξ(sol,0,1)[idm], Ξ(tm), atol=atol_u)
        ## anti-rotating dynamic phasor
        @test isapprox(TFT.ψ(sol,0,1)[idm], Ψ(tm), atol=atol_u)
        @test isapprox(TFT.ψ(sol,1,1)[idm], TFT.ξ(sol,1,1)[idm] * exp(-im * ω * tm), atol=atol_u)
        @test isapprox(TFT.ψ(sol,2,1)[idm], TFT.ξ(sol,2,1)[idm] * exp(-im * ω * tm), atol=atol_u)
        ## signal        
        @test isapprox(TFT.signal(sol)[idm], S(tm), atol=atol_u)
        @test isapprox(TFT.signal(sol,1)[idm], S(tm), atol=atol_u)
        ## error 
        @test isapprox(TFT.error(sol)[idm], 0.0u"V", atol=atol_u)
    end

end