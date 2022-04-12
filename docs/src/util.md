# Utilities

A number of functions are made available to the user to retrieve specific
components of the TFT solution

## Amplitude

```@docs
TFT.a(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)
```
```@docs
TFT.amplitude(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)
```

!!! note
    The amplitude `aₕ⁽¹⁾(t)` denotes the first derivate of the amplitude of the 
    zeroth derivative of the dynamic phasor `ξₕ⁽⁰⁾(t)`, not the amplitude of the 
    first derivative of the dynamic phasor `ξₕ⁽¹⁾(t)`.

## Phase

```@docs
TFT.ϕ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)
```
```@docs
TFT.phase(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)
```

!!! note
    The phase `ϕₕ⁽¹⁾(t)` denotes the first derivate of the phase of the 
    zeroth derivative of the dynamic phasor `ξₕ⁽⁰⁾(t)`, not the phase of the first 
    derivative of the dynamic phasor `ξₕ⁽¹⁾(t)`.

## Anti-Rotating Phase
```@docs
TFT.φ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)
```
```@docs
TFT.ar_phase(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)
```

!!! note
    The anti-rotating phase `φₕ⁽¹⁾(t)` denotes the first derivate of the 
    anti-rotating phase of the zeroth derivative of the dynamic phasor 
    `ξₕ⁽⁰⁾(t)`, not the anti-rotating phase of the first derivative of the 
    dynamic phasor `ξₕ⁽¹⁾(t)`.

## Frequency
```@docs
TFT.f(sol::TFT.AbstractDTFTSolution, H::Int=1)
```
```@docs
TFT.frequency(sol::TFT.AbstractDTFTSolution, H::Int=1)
```

## Dynamic Phasor
```@docs
TFT.ξ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)
```
```@docs
TFT.phasor(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)
```

## Anti-Rotating Dynamic Phasor
```@docs
TFT.ψ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)
```
```@docs
TFT.ar_phasor(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)
```

## Signal
```@docs
TFT.signal(sol::TFT.AbstractDTFTSolution)
```
```@docs
TFT.signal(sol::TFT.AbstractDTFTSolution, H::Int)
```

## Error
```@docs
TFT.error(sol::TFT.AbstractDTFTSolution)
```