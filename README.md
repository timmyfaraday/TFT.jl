# TFT.jl

<a href="https://github.com/timmyfaraday/TFT.jl/actions?query=workflow%3ACI"><img src="https://github.com/timmyfaraday/TFT.jl/workflows/CI/badge.svg"></img></a>
<a href="https://codecov.io/gh/timmyfaraday/TFT.jl"><img src="https://img.shields.io/codecov/c/github/timmyfaraday/TFT.jl?logo=Codecov"></img></a>
<a href="https://timmyfaraday.github.io/TFT.jl/"><img src="https://github.com/timmyfaraday/TFT.jl/workflows/Documentation/badge.svg"></img></a>

[![Active Development](https://img.shields.io/badge/Maintenance%20Level-Actively%20Developed-brightgreen.svg)](https://github.com/timmyfaraday/TFT.jl)

## Overview

TFT.jl is a Julia package for Taylor-Fourier transform. Similar to the Fourier transform, a Taylor-Fourier transform enables decomposing a time-variant function into its corresponding temporal frequency components. However, contrary to the Fourier transform, the Taylor-Fourier transform does not require the input signal to be periodic, rather approximates it aperiodicity through a Taylor polynomial. Consequently, the resulting frequency components are called dynamic phasors, rather than, as for the Fourier transform, simply phasors.

## Installation

The latest stable release of TFT.jl can be installed using the Julia package 
manager:

```julia
(v1.6) pkg> add "https://github.com/timmyfaraday/TFT.jl"
```
This package supports Julia v1.6 and later.

In order to test whether the package works, run:
```julia
(v1.6) pkg> test TFT
```
