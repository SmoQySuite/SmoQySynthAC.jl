@doc raw"""
    spectral_to_imaginary_time_correlation_function(;
        # KEYWORD ARGUMENTS
        τ::AbstractVector{T},
        β::T,
        spectral_function::Function,
        kernel_function::Function,
        tol::T = 1e-10,
    ) where {T<:AbstractFloat}

Calculate and return the imaginary-time correlation function
```math
C(\tau) = \int_{-\infty}^{\infty} d\omega \ K(\omega, \tau, \beta) \ A(\omega).
```
on a grid of ``\tau`` (`τ`) values, given a spectral function ``A(\omega)`` (`spectral_function`)
and kernel function ``K(\omega,\tau,\beta)`` (`kernel_function`). This integral is evaluated within
a specified tolerance `tol`.

## Arguments

- `τ::AbstractVector{T}`: Vector of imaginary time such that `τ[end] = β` equal the inverse temperature.
- `spectral_function::Function`: The spectral function ``A(\omega)`` that takes a single argument.
- `kernel_function::Function`: The kernel function ``K(\omega,\tau,\beta)`` that takes three arguments as shown.
- `tol::T = 1e-10`: Specified precision with which ``C(\tau)`` is evaluated.
"""
function spectral_to_imaginary_time_correlation_function(;
    # KEYWORD ARGUMENTS
    τ::AbstractVector{T},
    β::T,
    spectral_function::Function,
    kernel_function::Function,
    tol::T = 1e-10,
) where {T<:AbstractFloat}

    Cτ = quadgk(
        ω -> @.(kernel_function(ω, τ, β) * spectral_function(ω)),
        -Inf, 0.0, +Inf,
        atol = tol
    )[1]

    return Cτ
end


@doc raw"""
    spectral_to_matsubara_correlation_function(;
        # KEYWORD ARGUMENTS
        n::AbstractVector{Int},
        β::T,
        spectral_function::Function,
        kernel_function::Function,
        tol::T = 1e-10,
    ) where {T<:AbstractFloat}

Calculate and return the Matsubara correlation function
```math
C({\rm i} \omega_n) = \int_{-\infty}^{\infty} d\omega \ K_\beta(\omega, \omega_n) \ A(\omega).
```
for a vector of ``n \in \mathbb{Z}`` values, given a spectral function ``A(\omega)`` (`spectral_function`)
and kernel function ``K(\omega,\omega_n,\beta)`` (`kernel_function`). This integral is evaluated within
a specified tolerance `tol`.

Note that the kernel function should be called as `kernel_function(ω, n, β)` where `n`
specifies the Matsubara frequency, which is evaluated internally as either ``\omega_n = (2n+1)\pi/\beta``
or ``\omega_n = 2n\pi/\beta`` depending on whether the kernel function is fermionic of bosonic respectively.

## Arguments

- `n::AbstractVector{Int}`: Vector of integers specifying Matsubara frequencies for which ``C({\rm i}\omega_n) will be evaluated.
- `spectral_function::Function`: The spectral function ``A(\omega)`` that takes a single argument.
- `kernel_function::Function`: The kernel function ``K(\omega,\omega_n,\beta)`` that takes three arguments as shown.
- `tol::T = 1e-10`: Specified precision with which ``C({\rm i}\omega_n)`` is evaluated.
"""
function spectral_to_matsubara_correlation_function(;
    # KEYWORD ARGUMENTS
    n::AbstractVector{Int},
    β::T,
    spectral_function::Function,
    kernel_function::Function,
    tol::T = 1e-10,
) where {T<:AbstractFloat}

    Cτ = quadgk(
        ω -> @.(kernel_function(ω, n, β) * spectral_function(ω)),
        -Inf, 0.0, +Inf,
        atol = tol
    )[1]

    return Cτ
end