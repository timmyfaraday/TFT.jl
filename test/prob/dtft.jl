################################################################################
#  Copyright 2022, Tom Van Acker (BASF Antwerp)                                #
################################################################################
# TFT.jl                                                                       #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TFT.jl                                    #
################################################################################

@testset "DTFT" begin
    
    @testset "Fundamental Periodic Signal" begin 
        # amplitude 
        A(t)    = 10.0
        dA(t)   = 0.0
        d2A(t)  = 0.0
        
        # angle
        Φ(t)    = pi / 2
        dΦ(t)   = 0.0
        d2Φ(t)  = 0.0

        # analytic phasor, manually derived 
        P(t)    = A.(t) .* exp.(im .* Φ.(t))
        dP(t)   = (dA.(t) .+ im .* dΦ.(t) .* A.(t)) .* exp.(im .* Φ.(t))
        d2P(t)  = ((d2A.(t) .- dΦ.(t).^2 .* A.(t)) .+ im .* (2.0 .* dΦ.(t) .* dA.(t) .+ d2Φ.(t) .* A.(t))) .* exp.(im .* Φ.(t))            

        # analytic signal, manually derived
        S(t)    = real(P.(t) .* exp.(im .* ω .* t))
        dS(t)   = real((dP.(t) .+ im .* ω .* P.(t)) .* exp.(im .* ω .* t))
        d2S(t)  = real((d2P.(t) .- ω^2 .* P.(t)) .+ im .* 2.0 .* ω .* dP.(t) .* exp.(im .* ω .* t))

        # perform taylor-fourier transform
        sol = tft(S.(t), collect(t), [1], D, F, K)

        # amplitude tests
        @test isapprox(TFT.a(sol,0,1)[idm], A(tm), atol=atol)
        @test isapprox(TFT.a(sol,1,1)[idm], dA(tm), atol=atol)
        @test isapprox(TFT.a(sol,2,1)[idm], d2A(tm), atol=atol)
        # angle tests
        @test isapprox(TFT.φ(sol,0,1)[idm], Φ(tm), atol=atol)
        @test isapprox(TFT.φ(sol,1,1)[idm], dΦ(tm), atol=atol)
        @test isapprox(TFT.φ(sol,2,1)[idm], d2Φ(tm), atol=atol)
        # frequency test
        @test isapprox(TFT.frequency(sol,1)[idm], F + TFT.φ(sol,1,1)[idm] / (2 * pi), atol=atol)
        # rocof test
        @test isapprox(TFT.rocof(sol,1)[idm], TFT.φ(sol,2,1)[idm] / (2 * pi)^2, atol=atol)
        # phasor tests 
        @test isapprox(TFT.phasor(sol,0,1)[idm], P(tm), atol=atol)
        @test isapprox(TFT.phasor(sol,1,1)[idm], dP(tm), atol=atol)
        @test isapprox(TFT.phasor(sol,2,1)[idm], d2P(tm), atol=atol)
        # signal tests
        @test isapprox(TFT.signal(sol,0,1)[idm], S(tm), atol=atol)
        @test isapprox(TFT.signal(sol,1,1)[idm], dS(tm), atol=atol)
        @test isapprox(TFT.signal(sol,2,1)[idm], d2S(tm), atol=atol)
        # overall signal tests
        @test isapprox(TFT.overall_signal(sol,0)[idm], S(tm), atol=atol)
        @test isapprox(TFT.overall_signal(sol,1)[idm], dS(tm), atol=atol)
        @test isapprox(TFT.overall_signal(sol,2)[idm], d2S(tm), atol=atol)
        # error tests
        @test isapprox(TFT.error(sol)[idm], 0.0, atol=atol)
    end

    @testset "Fundamental Aperiodic Linear" begin 
        # amplitude 
        A(t)    = 10.0 - t
        dA(t)   = - 1.0
        d2A(t)  = 0.0

        # angle
        Φ(t)    = pi / 2 * t
        dΦ(t)   = pi / 2
        d2Φ(t)  = 0.0

        # analytic phasor, manually derived 
        P(t)    = A.(t) .* exp.(im .* Φ.(t))
        dP(t)   = (dA.(t) .+ im .* dΦ.(t) .* A.(t)) .* exp.(im .* Φ.(t))
        d2P(t)  = ((d2A.(t) .- dΦ.(t).^2 .* A.(t)) .+ im .* (2.0 .* dΦ.(t) .* dA.(t) .+ d2Φ.(t) .* A.(t))) .* exp.(im .* Φ.(t))            

        # analytic signal, manually derived
        S(t)    = real(P.(t) .* exp.(im .* ω .* t))
        dS(t)   = real((dP.(t) .+ im .* ω .* P.(t)) .* exp.(im .* ω .* t))
        d2S(t)  = real((d2P.(t) .- ω^2 .* P.(t)) .+ im .* 2.0 .* ω .* dP.(t) .* exp.(im .* ω .* t))

        # perform taylor-fourier transform
        sol = tft(S.(t), collect(t), [1], D, F, K)

        # amplitude tests
        @test isapprox(TFT.a(sol,0,1)[idm], A(tm), atol=atol)
        @test isapprox(TFT.a(sol,1,1)[idm], dA(tm), atol=atol)
        @test isapprox(TFT.a(sol,2,1)[idm], d2A(tm), atol=atol)
        # angle tests
        @test isapprox(TFT.φ(sol,0,1)[idm], Φ(tm), atol=atol)
        @test isapprox(TFT.φ(sol,1,1)[idm], dΦ(tm), atol=atol)
        @test isapprox(TFT.φ(sol,2,1)[idm], d2Φ(tm), atol=atol)
        # frequency test
        @test isapprox(TFT.frequency(sol,1)[idm], F + TFT.φ(sol,1,1)[idm] / (2 * pi), atol=atol)
        # rocof test
        @test isapprox(TFT.rocof(sol,1)[idm], TFT.φ(sol,2,1)[idm] / (2 * pi)^2, atol=atol)
        # phasor tests 
        @test isapprox(TFT.phasor(sol,0,1)[idm], P(tm), atol=atol)
        @test isapprox(TFT.phasor(sol,1,1)[idm], dP(tm), atol=atol)
        @test isapprox(TFT.phasor(sol,2,1)[idm], d2P(tm), atol=atol)
        # signal tests
        @test isapprox(TFT.signal(sol,0,1)[idm], S(tm), atol=atol)
        @test isapprox(TFT.signal(sol,1,1)[idm], dS(tm), atol=atol)
        @test isapprox(TFT.signal(sol,2,1)[idm], d2S(tm), atol=atol)
        # overall signal tests
        @test isapprox(TFT.overall_signal(sol,0)[idm], S(tm), atol=atol)
        @test isapprox(TFT.overall_signal(sol,1)[idm], dS(tm), atol=atol)
        @test isapprox(TFT.overall_signal(sol,2)[idm], d2S(tm), atol=atol)
        # error tests
        @test isapprox(TFT.error(sol)[idm], 0.0, atol=atol)
    end

    @testset "Fundamental Aperiodic Quadratic" begin 
        # amplitude
        A(t)    = 10.0 - t + 2.0 .* t^2
        dA(t)   = - 1.0 + 4.0 .* t
        d2A(t)  = 4.0

        # angle
        Φ(t)    = pi / 2 * t^2
        dΦ(t)   = pi * t
        d2Φ(t)  = pi

        # analytic phasor, manually derived 
        P(t)    = A.(t) .* exp.(im .* Φ.(t))
        dP(t)   = (dA.(t) .+ im .* dΦ.(t) .* A.(t)) .* exp.(im .* Φ.(t))
        d2P(t)  = ((d2A.(t) .- dΦ.(t).^2 .* A.(t)) .+ im .* (2.0 .* dΦ.(t) .* dA.(t) .+ d2Φ.(t) .* A.(t))) .* exp.(im .* Φ.(t))            

        # analytic signal, manually derived
        S(t)    = real(P.(t) .* exp.(im .* ω .* t))
        dS(t)   = real((dP.(t) .+ im .* ω .* P.(t)) .* exp.(im .* ω .* t))
        d2S(t)  = real((d2P.(t) .- ω^2 .* P.(t)) .+ im .* 2.0 .* ω .* dP.(t) .* exp.(im .* ω .* t))

        # perform taylor-fourier transform
        sol = tft(S.(t), collect(t), [1], D, F, K)

        # amplitude tests
        @test isapprox(TFT.a(sol,0,1)[idm], A(tm), atol=atol)
        @test isapprox(TFT.a(sol,1,1)[idm], dA(tm), atol=atol)
        @test isapprox(TFT.a(sol,2,1)[idm], d2A(tm), atol=atol)
        # angle tests
        @test isapprox(TFT.φ(sol,0,1)[idm], Φ(tm), atol=atol)
        @test isapprox(TFT.φ(sol,1,1)[idm], dΦ(tm), atol=atol)
        @test isapprox(TFT.φ(sol,2,1)[idm], d2Φ(tm), atol=atol)
        # frequency test
        @test isapprox(TFT.frequency(sol,1)[idm], F + TFT.φ(sol,1,1)[idm] / (2 * pi), atol=atol)
        # rocof test
        @test isapprox(TFT.rocof(sol,1)[idm], TFT.φ(sol,2,1)[idm] / (2 * pi)^2, atol=atol)
        # phasor tests 
        @test isapprox(TFT.phasor(sol,0,1)[idm], P(tm), atol=atol)
        @test isapprox(TFT.phasor(sol,1,1)[idm], dP(tm), atol=atol)
        @test isapprox(TFT.phasor(sol,2,1)[idm], d2P(tm), atol=atol)
        # signal tests
        @test isapprox(TFT.signal(sol,0,1)[idm], S(tm), atol=atol)
        @test isapprox(TFT.signal(sol,1,1)[idm], dS(tm), atol=atol)
        @test isapprox(TFT.signal(sol,2,1)[idm], d2S(tm), atol=atol)
        # overall signal tests
        @test isapprox(TFT.overall_signal(sol,0)[idm], S(tm), atol=atol)
        @test isapprox(TFT.overall_signal(sol,1)[idm], dS(tm), atol=atol)
        @test isapprox(TFT.overall_signal(sol,2)[idm], d2S(tm), atol=atol)
        # error tests
        @test isapprox(TFT.error(sol)[idm], 0.0, atol=atol)
    end

end