########################################
## GENERIC MATSUBARA KERNEL  FUNCTION ##
########################################


@doc raw"""
    kernel_mat(ω::T, ω_n::T) where {T<:AbstractFloat}

Generic Matsubara frequency kernel
```math
K_\beta(\omega, \omega_n) = \frac{1}{\omega - {\rm i}\omega_n},
```
where ``\omega_n = (2n+1)\pi/\beta \ (= 2n\pi/\beta)`` for fermions (bosons) with ``n \in \mathbb{Z}``.
This kernel assumes the sign convention
```math
C(\tau) = \langle \hat{a}(\tau) \hat{a}^\dagger(0) \rangle
```
for a the Green's function, where it is assumed that ``\tau \in [0,\beta)``.
"""
kernel_mat(ω::T, ω_n::T) where {T<:AbstractFloat} = inv( ω - im * ω_n )


###############################
## FERMIONIC KERNEL FUNCTION ##
###############################


@doc raw"""
    kernel_tau_fermi(ω::T, τ::T, β::T) where {T<:AbstractFloat}

The imaginary time fermionic kernel
```math
\begin{align*}
K_\beta(\omega,\tau) & = \overbrace{\left(\frac{e^{-\tau\omega}}{1+e^{-\beta\omega}}\right)}^{\text{numerically unstable}} \\
                     & = \underbrace{\left( e^{\tau\omega} + e^{(\tau-\beta)\omega} \right)^{-1}}_{\text{numerically stable}},
\end{align*}
```
where it is assumed that ``\tau \in [0,\beta)``.
"""
kernel_tau_fermi(ω::T, τ::T, β::T) where {T<:AbstractFloat} = inv(exp(τ*ω) + exp((τ-β)*ω))


@doc raw"""
    kernel_mat_fermi(ω::T, n::Int, β::T) where {T<:AbstractFloat}

The fermionic matsubara frequency kernel
```math
K_\beta(\omega, \omega_n) = \frac{1}{{\omega - {\rm i}\omega_n}},
```
where ``\omega_n = (2n+1)\pi/\beta`` for fermions with ``n \in \mathbb{Z}``.
"""
kernel_mat_fermi(ω::T, n::Int, β::T) where {T<:AbstractFloat} = kernel_mat(ω, (2*n+1)*π/β)


##############################
## BOSONIC KERNEL FUNCTIONS ##
##############################


@doc raw"""
    kernel_tau_bose(ω::T, τ::T, β::T) where {T<:AbstractFloat}

The imaginary time bosonic kernel
```math
\begin{align*}
K_\beta(\omega,\tau) & = \overbrace{\left(\frac{e^{-\tau\omega}}{1-e^{-\beta\omega}}\right)}^{\text{numerically unstable}} \\
                     & = \underbrace{\left( e^{\tau\omega} - e^{(\tau-\beta)\omega} \right)^{-1}}_{\text{numerically stable}},
\end{align*}
```
where it is assumed that ``\tau \in [0,\beta)``.
"""
kernel_tau_bose(ω::T, τ::T, β::T) where {T<:AbstractFloat} = inv(exp(τ*ω) - exp((τ-β)*ω))


@doc raw"""
    kernel_mat_bose(ω::T, n::Int, β::T) where {T<:AbstractFloat}

The bosonic matsubara frequency kernel
```math
K_\beta(\omega, \omega_n) = \frac{1}{\omega - {\rm i}\omega_n},
```
where ``\omega_n = 2n\pi/\beta`` for bosons with ``n \in \mathbb{Z}``.
"""
kernel_mat_bose(ω::T, n::Int, β::T) where {T<:AbstractFloat} = kernel_mat(ω, 2*n*π/β)


#######################################
## MODIFIED BOSONIC KERNEL FUNCTIONS ##
#######################################


@doc raw"""
    kernel_tau_bose_alt(ω::T, τ::T, β::T) where {T<:AbstractFloat}

The imaginary time bosonic kernel
```math
\begin{align*}
K_\beta(\omega,\tau) & = \overbrace{\left(\frac{\omega e^{-\tau\omega}}{1-e^{-\beta\omega}}\right)}^{\text{numerically unstable}} \\
                     & = \underbrace{\omega \left( e^{\tau\omega} - e^{(\tau-\beta)\omega} \right)^{-1}}_{\text{numerically stable}},
\end{align*}
```
where it is assumed that ``\tau \in [0,\beta)`` and ``K_\beta(0,\tau) = \beta^{-1}``.
"""
kernel_tau_bose_alt(ω::T, τ::T, β::T) where {T<:AbstractFloat} = iszero(ω) ? inv(β) : ω * kernel_tau_bose(ω, τ, β)


@doc raw"""
    kernel_mat_bose_alt(ω::T, n::Int, β::T) where {T<:AbstractFloat}

The bosonic matsubara frequency kernel
```math
K_\beta(\omega, \omega_n) = \frac{\omega}{\omega - {\rm i}\omega_n},
```
where ``\omega_n = 2n\pi/\beta`` for bosons with ``n \in \mathbb{Z}`` and ``K_\beta(0,0) = -1``.
"""
kernel_mat_bose_alt(ω::T, n::Int, β::T) where {T<:AbstractFloat} = (iszero(ω) && iszero(n)) ? -one(T) : ω * kernel_mat_bose(ω, n, β)


########################################
## SYMMETRIC BOSONIC KERNEL FUNCTIONS ##
########################################


@doc raw"""
    kernel_tau_sym_bose(ω::T, τ::T, β::T) where {T<:AbstractFloat}

The imaginary time symmetrized bosonic kernel
```math
\begin{align*}
K_\beta(\omega,\tau) & = \overbrace{\left(\frac{e^{-\tau\omega} + e^{-(\beta-\tau)\omega}}{1-e^{-\beta\omega}}\right)}^\text{numerically unstable} \\
& = \underbrace{\left( e^{\tau\omega} - e^{(\tau-\beta)\omega} \right)^{-1} - \left(e^{-\tau\omega} - e^{-(\tau-\beta)\omega} \right)^{-1}}_{\text{numerically stable}},
\end{align*}
```
where it is assumed that ``\tau \in [0,\beta)``.
"""
kernel_tau_sym_bose(ω::T, τ::T, β::T) where {T<:AbstractFloat} = inv(exp(τ*ω) - exp((τ-β)*ω)) - inv(exp(-τ*ω) - exp(-(τ-β)*ω))


@doc raw"""
    kernel_mat_sym_bose(ω::T, n::Int, β::T) where {T<:AbstractFloat}

The symmetrized bosonic matsubara frequency kernel
```math
\begin{align*}
K_β(\omega, \omega_n) & = \frac{2\omega}{\omega_n^2 + \omega^2},
\end{align*}
```
where ``\omega_n = 2n\pi/\beta`` for bosons with ``n \in \mathbb{Z}``.
"""
kernel_mat_sym_bose(ω::T, n::Int, β::T) where {T<:AbstractFloat} = kernel_mat_sym_bose(ω, 2*n*π/β)
kernel_mat_sym_bose(ω::T, ω_n::T) where {T<:AbstractFloat} = 2*ω/(ω_n^2 + ω^2)


#################################################
## MODIFIED SYMMETRIC BOSONIC KERNEL FUNCTIONS ##
#################################################


@doc raw"""
    kernel_tau_sym_bose_alt(ω::T, τ::T, β::T) where {T<:AbstractFloat}

The imaginary time symmetrized bosonic kernel
```math
\begin{align*}
K_\beta(\omega,\tau) & = \overbrace{\omega \left(\frac{e^{-\tau\omega} + e^{-(\beta-\tau)\omega}}{1-e^{-\beta\omega}}\right)}^\text{numerically unstable} \\
& = \underbrace{\omega\left[\left( e^{\tau\omega} - e^{(\tau-\beta)\omega} \right)^{-1} - \left(e^{-\tau\omega} - e^{-(\tau-\beta)\omega} \right)^{-1}\right]}_{\text{numerically stable}},
\end{align*}
```
where it is assumed that ``\tau \in [0,\beta)`` and ``K_\beta(0,\tau) = 2/\beta``.
"""
kernel_tau_sym_bose_alt(ω::T, τ::T, β::T) where {T<:AbstractFloat} = izero(ω) ? 2/β : ω * kernel_tau_sym_bose(ω, τ, β)


@doc raw"""
    kernel_mat_sym_bose_alt(ω::T, n::Int, β::T) where {T<:AbstractFloat}

The symmetrized bosonic matsubara frequency kernel
```math
\begin{align}
K_β(\omega, \omega_n) & = \frac{2\omega^2}{\omega_n^2 + \omega^2},
\end{align}
```
where ``\omega_n = 2n\pi/\beta`` for bosons with ``n \in \mathbb{Z}`` and ``K_\beta(0,0) = 2``.
"""
kernel_mat_sym_bose_alt(ω::T, n::Int, β::T) where {T<:AbstractFloat} = (iszero(ω) && iszero(n)) ? 2.0 : ω * kernel_mat_sym_bose(ω, n, β)


##################################
## QUANTUM STATISTICS FUNCTIONS ##
##################################


@doc raw"""
    fermi(ϵ::T, β::T) where {T<:AbstractFloat}

The Fermi-Dirac function
```math
\begin{align*}
f_\beta(\epsilon) & = \overbrace{\left( \frac{1}{e^{\beta\epsilon} + 1}\right)}^{\text{numerically unstable}} \\
                  & = \underbrace{\frac{1}{2}\left( 1 - \tanh\left(\frac{\beta\epsilon}{2}\right) \right)}_{\text{numerically stable}},
\end{align*}
```
where ``\epsilon`` is energy and ``\beta`` is inverse temperature.
"""
fermi(ϵ::T, β::T) where {T<:AbstractFloat} = (1 - tanh(β*ϵ/2))/2


@doc raw"""
    bose(ϵ::T, β::T) where {T<:AbstractFloat}

The Bose-Einstein function
```math
\begin{align*}
n_\beta(\epsilon) & = \overbrace{\left( \frac{1}{e^{\beta\epsilon} - 1}\right)}^{\text{numerically unstable}} \\
                  & = \underbrace{\frac{1}{2}\left(\coth\left(\frac{\beta\epsilon}{2}\right) - 1 \right)}_{\text{numerically stable}},
\end{align*}
```
where ``\epsilon`` is energy and ``\beta`` is inverse temperature.
"""
bose(ϵ::T, β::T) where {T<:AbstractFloat} = (coth(β*ϵ/2) - 1)/2