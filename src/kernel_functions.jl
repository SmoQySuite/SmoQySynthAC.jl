@doc raw"""
    kernel_tau_fermi(ω::T, τ::T, β::T) where {T<:AbstractFloat}

The imaginary time fermionic kernel
```math
\begin{align}
K_\beta(\omega,\tau) & = \overbrace{\left(\frac{e^{-\tau\omega}}{1+e^{-\beta\omega}}\right)}^{\text{numerically unstable}} \\
                     & = \underbrace{\left( e^{\tau\omega} + e^{(\tau-\beta)\omega} \right)^{-1}}_{\text{numerically stable}},
\end{align}
```
where it is assumed that ``\tau \in [0,\beta)``.
"""
kernel_tau_fermi(ω::T, τ::T, β::T) where {T<:AbstractFloat} = inv(exp(τ*ω) + exp((τ-β)*ω))
# kernel_tau_fermi(ω::T, τ::T, β::T) where {T<:AbstractFloat} = (-β*ω > 100) ? exp((β-τ)*ω) : exp(-τ*ω)/(1.0 + exp(-β*ω))


@doc raw"""
    kernel_tau_bose(ω::T, τ::T, β::T) where {T<:AbstractFloat}

The imaginary time bosonic kernel
```math
\begin{align}
K_\beta(\omega,\tau) & = \overbrace{\left(\frac{e^{-\tau\omega}}{1-e^{-\beta\omega}}\right)}^{\text{numerically unstable}} \\
                     & = \underbrace{\left( e^{\tau\omega} - e^{(\tau-\beta)\omega} \right)^{-1}}_{\text{numerically stable}},
\end{align}
```
where it is assumed that ``\tau \in [0,\beta)``.
"""
kernel_tau_bose(ω::T, τ::T, β::T) where {T<:AbstractFloat} = inv(exp(τ*ω) - exp((τ-β)*ω))
# kernel_tau_bose(ω::T, τ::T, β::T) where {T<:AbstractFloat} = (-β*ω > 100) ? -exp((β-τ)*ω) : -exp(-τ*ω)/expm1(-β*ω)


@doc raw"""
    kernel_tau_sym_bose(ω::T, τ::T, β::T) where {T<:AbstractFloat}

The imaginary time symmetrized bosonic kernel
```math
\begin{align}
K_\beta(\omega,\tau) & = \overbrace{\left(\frac{e^{-\tau\omega} + e^{-(\beta-\tau)\omega}}{1-e^{-\beta\omega}}\right)}^\text{numerically unstable} \\
& = \underbrace{\left( e^{\tau\omega} - e^{(\tau-\beta)\omega} \right)^{-1} - \left(e^{-\tau\omega} - e^{-(\tau-\beta)\omega} \right)^{-1}}_{\text{numerically stable}},
\end{align}
```
where it is assumed that ``\tau \in [0,\beta)``.
"""
kernel_tau_sym_bose(ω::T, τ::T, β::T) where {T<:AbstractFloat} = inv(exp(τ*ω) - exp((τ-β)*ω)) - inv(exp(-τ*ω) - exp(-(τ-β)*ω))
# kernel_tau_sym_bose(ω::T, τ::T, β::T) where {T<:AbstractFloat} = (-β*ω > 100) ? -(exp((β-τ)*ω) + exp(τ*ω)) : -(exp(-τ*ω) + exp(-(β-τ)*ω))/expm1(-β*ω)


@doc raw"""
    kernel_mat_fermi(ω::T, n::Int, β::T) where {T<:AbstractFloat}

The fermionic matsubara frequency kernel
```math
K_\beta(\omega, \omega_n) = \frac{1}{{\rm i}\omega_n - \omega},
```
where ``\omega_n = (2n+1)\pi/\beta`` for fermions with ``n \in \mathbb{Z}``.
"""
kernel_mat_fermi(ω::T, n::Int, β::T) where {T<:AbstractFloat}  = kernel_mat(ω, (2*n+1)*π/β)


@doc raw"""
    kernel_mat_bose(ω::T, n::Int, β::T) where {T<:AbstractFloat}

The bosonic matsubara frequency kernel
```math
K_\beta(\omega, \omega_n) = \frac{1}{{\rm i}\omega_n - \omega},
```
where ``\omega_n = 2n\pi/\beta`` for bosons with ``n \in \mathbb{Z}``.
"""
kernel_mat_bose(ω::T, n::Int, β::T) where {T<:AbstractFloat}  = kernel_mat(ω, 2*n*π/β)

# matsubara kernel function
kernel_mat(ω::T, ω_n::T) where {T<:AbstractFloat} = inv( im * ω_n - ω )

@doc raw"""
    fermi(ϵ::T, β::T) where {T<:AbstractFloat}

The Fermi-Dirac function
```math
\begin{align}
f_\beta(\epsilon) & = \overbrace{\left( \frac{1}{e^{\beta\epsilon} + 1}\right)}^{\text{numerically unstable}} \\
                  & = \underbrace{\frac{1}{2}\left( 1 - \tanh\left(\frac{\beta\epsilon}{2}\right) \right)}_{\text{numerically stable}},
\end{align}
```
where ``\epsilon`` is energy and ``\beta`` is inverse temperature.
"""
fermi(ϵ::T, β::T) where {T<:AbstractFloat} = (1 - tanh(β*ϵ/2))/2

@doc raw"""
    bose(ϵ::T, β::T) where {T<:AbstractFloat}

The Bose-Einstein function
```math
\begin{align}
n_\beta(\epsilon) & = \overbrace{\left( \frac{1}{e^{\beta\epsilon} - 1}\right)}^{\text{numerically unstable}} \\
                  & = \underbrace{\frac{1}{2}\left(\coth\left(\frac{\beta\epsilon}{2}\right) - 1 \right)}_{\text{numerically stable}},
\end{align}
```
where ``\epsilon`` is energy and ``\beta`` is inverse temperature.
"""
bose(ϵ::T, β::T) where {T<:AbstractFloat} = (coth(β*ϵ/2) - 1)/2