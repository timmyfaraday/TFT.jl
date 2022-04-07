# Utilities

A number of functions are made available to the user to retrieve specific
components of the TFT solution

## Quick Links

- [Utilities](#utilities)
  - [Quick Links](#quick-links)
  - [Amplitude](#amplitude)
  - [Phase](#phase)
  - [Anti-Rotating Phase](#anti-rotating-phase)
  - [Dynamic Phasor](#dynamic-phasor)
  - [Anti-Rotating Dynamic Phasor](#anti-rotating-dynamic-phasor)
  - [Signal](#signal)
  - [Error](#error)

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