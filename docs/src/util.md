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

## Phase

```@docs
TFT.ϕ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)
```
```@docs
TFT.phase(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)
```

## Anti-Rotating Phase
```@docs
TFT.φ(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)
```
```@docs
TFT.ar_phase(sol::TFT.AbstractDTFTSolution, D::Int=0, H::Int=1)
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