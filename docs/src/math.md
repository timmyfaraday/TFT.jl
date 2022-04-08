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

### Zeroth Derivative

The zeroth derivative of the hth-harmonic dynamic phasor `ξₕ⁽⁰⁾(t)` is given by,
```math
\begin{aligned}
    \xi^{(0)}_{h}(t)        =& \frac{a^{(0)}_{h}(t)}{2} \, \exp(j \phi^{(0)}_{h}(t)) \\
                            =& \frac{a^{(0)}_{h}(t)}{2} \, \exp(j \varphi^{(0)}_{h}(t)) \, \exp(j \omega_{h} t) \\
                            =& \psi^{(0)}_{h}(t) \, \exp(j \omega_{h} t)

\end{aligned}
```
where `ψₕ⁽⁰⁾(t)` denotes the anti-rotating zeroth derivative of the hth-harmonic 
dynamic phasor. 

The zeroth derivative of the amplitude `aₕ⁽⁰⁾(t)`, phase `ϕₕ⁽⁰⁾(t)` and anti-rotating 
phase `φₕ⁽⁰⁾(t)` are, respectively,
```math
\begin{aligned}
    a^{(0)}_{h}(t)          =& |2 \, \xi^{(0)}_{h}(t)| \\
    \phi^{(0)}_{h}(t)       =& \angle \big( \xi^{(0)}_{h}(t) \big) \\
    \varphi^{(0)}_{h}(t)    =& \angle \big( \xi^{(0)}_{h}(t) \, \exp(-j \omega_{h} t) \big) = \angle \big( \psi^{(0)}_{h}(t) \big) = \phi^{(0)}_{h}(t) - \omega_{h}t
\end{aligned}
```
The zeroth derivative hth-harmonic signal `sₕ⁽⁰⁾(t)` and overall signal `s⁽⁰⁾(t)` 
are, respectively,
```math
\begin{aligned}
    s^{(0)}_{h}(t)          =& \xi^{(0)}_{h}(t) + \text{conj}\big( \xi^{(0)}_{h}(t) \big) \\
    s^{(0)}(t)              =& \sum_{h \in H} \xi^{(0)}_{h}(t) + \text{conj}\big( \xi^{(0)}_{h}(t) \big),
\end{aligned}
```
where `conj(ξₕ⁽⁰⁾(t))` denotes the complex conjugate of the zeroth derivative of the 
hth-harmonic dynamic phasor.

### First Derivative

The first derivate of the hth-harmonic dynamic phasor `ξₕ⁽¹⁾(t)` is given by,
```math
\begin{aligned}
    \xi_{h}^{(1)}(t)        =&  \frac{\mathrm{d}\xi_{h}^{(0)}(t)}{\mathrm{d}t} \\
                            =&  \frac{1}{2} \big( a^{(1)}_{h}(t) +
                                j \, \phi^{(1)}_{h}(t) \, a^{(0)}_{h}(t)) \big)
                                \, \exp(j \phi^{(0)}_{h}(t))

\end{aligned}
```

The first derivative of the hth-harmomic anti-rotating dynamic phasor `ψₕ⁽¹⁾(t)`
is,
```math 
\begin{aligned}
    \psi_{h}^{(1)}(t)       =&  \frac{\mathrm{d}\psi_{h}^{(0)}(t)}{\mathrm{d}t} \\
                            =&  \frac{\mathrm{d}\big(\xi_{h}^{(0)}(t) 
                                \exp(-j \omega_{h} t)\big)}{\mathrm{d}t} \\
                            =&  \big(\xi_{h}^{(1)}(t) - j \omega_{h}
                                \xi_{h}^{(0)}(t) \big) \exp(-j \omega_{h} t)
\end{aligned}
```

The first derivative of the hth-harmonic amplitude `aₕ⁽¹⁾(t)`, phase `ϕₕ⁽¹⁾(t)` and
anti-rotating phase `φₕ⁽¹⁾(t)` are, respectively,
```math
\begin{aligned}
    a_{h}^{(1)}(t)          =& ℜ[2 \, \xi_{h}^{(1)}(t) \, \exp(-j \phi_{h}^{(0)}(t))] \\
    \phi_{h}^{(1)}(t)       =& \frac{ℑ[2 \, \xi_{h}^{(1)}(t) \, \exp(-j \phi_{h}^{(0)}(t))]}{a^{(0)}_{h}(t)} \\
    \varphi_{h}^{(1)}(t)    =& \phi_{h}^{(1)}(t) - \omega_{h}
\end{aligned}
```

The frequency `fₕ(t)` of the zeroth derivative hth-harmonic signal `sₕ⁽⁰⁾(t)` is,
```math
\begin{aligned}
    f_{h}(t)                =& TODO
\end{aligned}
```
