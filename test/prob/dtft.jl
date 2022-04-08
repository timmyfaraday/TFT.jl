################################################################################
#  Copyright 2022, Tom Van Acker (BASF Antwerp)                                #
################################################################################
# TFT.jl                                                                       #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TFT.jl                                    #
################################################################################

@testset "DTFT" begin
    
    @testset "Fundamental Periodic Signal" begin 
        # input
        A(t)    = 10.0 
        Φ(t)    = pi / 2

        # derived input
        Ξ(t)    = A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t)
        Ψ(t)    = A.(t) ./ 2 .* exp.(im .* Φ.(t))
        S(t)    = real.(A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t) .+ 
                  conj.(A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t)))

        # perform taylor-fourier transform
        sol = tft(S.(t), collect(t), [1], D, F, K)

        # tests
        @test isapprox(TFT.a(sol,0,1)[idm], A(tm), atol=atol)
        @test isapprox(TFT.a(sol,1,1)[idm], _FD.gradient(A,tm), atol=atol)
        @test isapprox(TFT.φ(sol,0,1)[idm], Φ(tm), atol=atol)
        @test isapprox(TFT.ξ(sol,0,1)[idm], Ξ(tm), atol=atol)
        @test isapprox(TFT.ψ(sol,0,1)[idm], Ψ(tm), atol=atol)
        @test isapprox(TFT.signal(sol)[idm], S(tm), atol=atol)
        @test isapprox(TFT.signal(sol,1)[idm], S(tm), atol=atol)
        @test isapprox(TFT.error(sol)[idm], 0.0, atol=atol)
    end

    @testset "Fundamental Aperiodic Linear" begin 
        # input
        A(t)    = 10.0 - t
        Φ(t)    = pi / 2 * t

        # derived input
        Ξ(t)    = A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t)
        Ψ(t)    = A.(t) ./ 2 .* exp.(im .* Φ.(t))
        S(t)    = real.(A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t) .+ 
                  conj.(A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t)))

        # perform taylor-fourier transform
        sol = tft(S.(t), collect(t), [1], D, F, K)

        # tests
        @test isapprox(TFT.a(sol,0,1)[idm], A(tm), atol=atol)
        @test isapprox(TFT.φ(sol,0,1)[idm], Φ(tm), atol=atol)
        @test isapprox(TFT.ξ(sol,0,1)[idm], Ξ(tm), atol=atol)
        @test isapprox(TFT.ψ(sol,0,1)[idm], Ψ(tm), atol=atol)
        @test isapprox(TFT.signal(sol)[idm], S(tm), atol=atol)
        @test isapprox(TFT.signal(sol,1)[idm], S(tm), atol=atol)
        @test isapprox(TFT.error(sol)[idm], 0.0, atol=atol)
    end

    @testset "Fundamental Aperiodic Quadratic" begin 
        # input
        A(t)    = 10.0 - t + 2.0 .* t^2
        Φ(t)    = pi / 2 * t^2

        # derived input
        Ξ(t)    = A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t)
        Ψ(t)    = A.(t) ./ 2 .* exp.(im .* Φ.(t))
        S(t)    = real.(A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t) .+ 
                  conj.(A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t)))

        # perform taylor-fourier transform
        sol = tft(S.(t), collect(t), [1], D, F, K)

        # tests
        @test isapprox(TFT.a(sol,0,1)[idm], A(tm), atol=atol)
        @test isapprox(TFT.φ(sol,0,1)[idm], Φ(tm), atol=atol)
        @test isapprox(TFT.ξ(sol,0,1)[idm], Ξ(tm), atol=atol)
        @test isapprox(TFT.ψ(sol,0,1)[idm], Ψ(tm), atol=atol)
        @test isapprox(TFT.signal(sol)[idm], S(tm), atol=atol)
        @test isapprox(TFT.signal(sol,1)[idm], S(tm), atol=atol)
        @test isapprox(TFT.error(sol)[idm], 0.0, atol=atol)
    end

end