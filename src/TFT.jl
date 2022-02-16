module TFT

    # import pkg
    import DSP
    import Polynomials

    # pkg constants 
    const _DSP = DSP
    const _POL = Polynomials

    # include
    include("ospline.jl")
    include("estimator.jl")

    # export
    export sample_ospline
    export harmonic_state_estimator
    
end
