@doc raw"""
    kernel_tau_fermi(ω::T, τ::T, β::T) where {T<:AbstractFloat}

The imaginary time fermionic kernel
```math
K_\beta(\omega,\tau) = \frac{e^{-\tau\omega}}{1+e^{-\beta\omega}},
```
where it is assumed that ``\tau \in [0,\beta)``.
"""
kernel_tau_fermi(ω::T, τ::T, β::T) where {T<:AbstractFloat} = (-β*ω > 100) ? exp((β-τ)*ω) : exp(-τ*ω)/(1.0 + exp(-β*ω))


@doc raw"""
    kernel_tau_bose(ω::T, τ::T, β::T) where {T<:AbstractFloat}

The imaginary time bosonic kernel
```math
K_\beta(\omega,\tau) = \frac{e^{-\tau\omega}}{1-e^{-\beta\omega}},
```
where it is assumed that ``\tau \in [0,\beta)``.
"""
kernel_tau_bose(ω::T, τ::T, β::T) where {T<:AbstractFloat} = (-β*ω > 100) ? -exp((β-τ)*ω) : -exp(-τ*ω)/expm1(-β*ω)


@doc raw"""
    kernel_tau_sym_bose(ω::T, τ::T, β::T) where {T<:AbstractFloat}

The imaginary time symmetrized bosonic kernel
```math
K_\beta(\omega,\tau) = \frac{e^{-\tau\omega} + e^{-(\beta-\tau)\omega}}{1-e^{-\beta\omega}},
```
where it is assumed that ``\tau \in [0,\beta)``.
"""
kernel_tau_sym_bose(ω::T, τ::T, β::T) where {T<:AbstractFloat} = (-β*ω > 100) ? -(exp((β-τ)*ω) + exp(τ*ω)) : -(exp(-τ*ω) + exp(-(β-τ)*ω))/expm1(-β*ω)


@doc raw"""
    kernel_mat_fermi(ω::T, n::Int, β::T) where {T<:AbstractFloat}

The fermionic matsubara frequency kernel
```math
K_\beta(\omega, {\rm i}\omega_n) = \frac{1}{{\rm i}\omega_n - \omega},
```
where ``\omega_n = (2n+1)\pi/\beta`` for fermions with ``n \in \mathbb{Z}``.
"""
kernel_mat_fermi(ω::T, n::Int, β::T) where {T<:AbstractFloat}  = kernel_mat(ω, (2*n+1)*π/β)


@doc raw"""
    kernel_mat_bose(ω::T, n::Int, β::T) where {T<:AbstractFloat}

The bosonic matsubara frequency kernel
```math
K_\beta(\omega, {\rm i}\omega_n) = \frac{1}{{\rm i}\omega_n - \omega},
```
where ``\omega_n = 2n\pi/\beta`` for bosons with ``n \in \mathbb{Z}``.
"""
kernel_mat_bose(ω::T, n::Int, β::T) where {T<:AbstractFloat}  = kernel_mat(ω, 2*n*π/β)

# matsubara kernel function
kernel_mat(ω::T, ω_n::T) where {T<:AbstractFloat} = inv( im * ω_n - ω )