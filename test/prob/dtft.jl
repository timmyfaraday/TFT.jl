################################################################################
#  Copyright 2022, Tom Van Acker (BASF Antwerp)                                #
################################################################################
# TFT.jl                                                                       #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TFT.jl                                    #
################################################################################

@testset "DTFT" begin
    
    @testset "Fundamental " begin 
        A, Φ = 10
        S(t) = A .* exp.(im .* Φ) .* exp.(im .* ω .* t)
        
        sol = tft(S(t), t, [1], 2)

        @test amplitude(sol,0) ≈ A
        @test amplitude(sol,1) ≈ 
        @test amplitude(sol,2) ≈ 

        @test angle(sol,0) 
    end

end