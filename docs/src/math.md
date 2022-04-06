# Mathematical Background on Taylor-Fourier Transform

## Nomenclature

| Symbol    | Description                   |
|:----------|:------------------------------|
| `s`       | signal                        |
| `a`       | amplitude                     |
| `ϕ`       | phase [(rad)]                 |
| `φ`       | anti-rotating phase [(rad)]   |
| `ξ`       | dynamic phasor                |
| `ψ`       | anti-rotating dynamic phasor  | 

| Set       | Description                   |
|:----------|:------------------------------|
| `h ∈ H`   | set of harmonic number [Int]  |

| Parameter | Description                   |
|:----------|:------------------------------|
| `F`       | frequency [Hz]                |
| `ω`       | angular frequency [(rad)/s]   |

## Taylor-Fourier Transform

The Taylor-Fourier Transform function `tft()` gives the up-to-Dth derivative of 
the hth-harmonic dynamic phasors `ξₕ⁽ᵈ⁾(t), ∀ d ∈ {0,..,D}, h ∈ H`. 

The zeroth derivative of the hth-harmonic dynamic phasor `ξₕ⁽⁰⁾(t)` is given by,
```math
\begin{aligned}
    \xi^{(0)}_{h}(t)        =& \frac{a^{(0)}_{h}(t)}{2} \, \exp(j \phi^{(0)}_{h}(t))
                            =& \frac{a^{(0)}_{h}(t)}{2} \, \exp(j \varphi^{(0)}_{h}(t)) \, \exp(j \omega_{h} t)
                            =& \psi^{(0)}_{h}(t) \, \exp(j \omega_{h} t)

\end{aligned}
```
where `ψₕ⁽⁰⁾(t)` denotes the anti-rotating zeroth derivative of the hth-harmonic 
dynamic phasor. The corresponding amplitude `aₕ⁽⁰⁾(t)`, phase `ϕₕ⁽⁰⁾(t)` and anti-
rotating phase `φₕ⁽⁰⁾(t)` are, respectively,
```
\begin{aligned}
    a^{(0)})_{h}(t)         =& 2 \cdot |\xi^{(0)}_{h}(t)|
    \phi^{(0)}_{h}(t)       =& \angle \xi^{(0)}_{h}(t)
    \varphi^{(0)}_{h}(t)    =& \angle \xi^{(0)}_{h}(t) \, \exp(-j \omega_{h} t)
\end{aligned}
```
The corresponding hth-harmonic signal `sₕ⁽⁰⁾(t)` and overall signal `s⁽⁰⁾(t)` are, 
respectively,
```math
\begin{aligned}
    s^{(0)}_{h}(t)          =& \xi^{(0)}_{h}(t) + \bar{\xi}^{(0)}_{h}(t)^{*}
    s^{(0)}(t)              =& \sum_{h \in H} \xi^{(0)}_{h}(t) + \bar{\xi}^{(0)}_{h}(t)^{*},
\end{aligned}
```
where `ξₕ⁽⁰⁾(t)*` denotes the complex conjugate of the zeroth derivative of the 
hth-harmonic dynamic phasor.
