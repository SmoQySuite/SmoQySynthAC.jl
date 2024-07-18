@doc raw"""
    add_noise(;
        # KEYWORD ARGUMENTS
        Cτ_exact::AbstractVector{T},
        τ::AbstractVector{T},
        σ::T,
        ξ::T,
        sum_rule::Function = C0 -> 1 - C0,
        noise_type::Symbol = :TruncatedNormal
    ) where {T<:AbstractFloat}

    add_noise(
        # ARGUMENTS
        N_samples::Int;
        # KEYWORD ARGUMENTS
        Cτ_exact::AbstractVector{T},
        τ::AbstractVector{T},
        σ::T,
        ξ::T,
        sum_rule::Function = C0 -> 1 - C0,
        noise_type::Symbol = :TruncatedNormal
    ) where {T<:AbstractFloat}

Add noise to an imaginary time correlation function ``C(\tau)`` that is exponentially correlated in imaginary time.

## Arguments

- `N_samples` (optional): Number of random samples to generate. If passed this function returns a matrix where the columsn correspond to the different samples.

## Keyword Arguments

- `Cτ_exact::AbstractVector{T}`: Vector containing the exact values for ``C(\tau)``.
- `τ::AbstractVector{T}`: Vector specifying the imaginary time ``\tau`` grid that ``C(\tau)`` is evaluated on. Assumes that the last element equals the inverse temperature, i.e. `τ[end] = β`.
- `σ::T`: Standard deviation of the noise; controls the typical amplitude of the error.
- `ξ::T`: Correlation length associated with the noise in imaginary time.
- `sum_rule::Function = C0 -> 1 - C0`: Enforces sum rule, or bounday condition in imaginay time. Default behavior assumes a fermionic correlation function, enforcing that $C(\beta) = 1 - C(0)$.
- `noise_type::Symbol = :TruncatedNormal`: Distribution that the noise is sampled from prior to correlations being introduced in imaginary time. Available options are `noise_type ∈ (:TruncatedNormal, :Gamma, :Normal)`.

## Additional Information

If the correlation function ``C(\tau_i)`` the corresponding noisy correlation function is given by
```math
C_{\rm noisy}(\tau_i) = C(\tau_i) + \frac{\sum_j e^{-|\tau_j-\tau_i|/\xi} R_j}{\sum_j e^{-2|\tau_j-\tau_i|/\xi}},
```
where the sum is performed assuming periodic boundary conditions. If `noise_type = :Normal` then the random numbers of sampled
according to ``R_j \sim {\rm Normal}(0,\sigma)``. Otherwise,
```math
R_j \sim {\rm sign}(C(\tau_j)) \cdot P(|C(\tau_j)|, \sigma) - C(\tau_j)
```
where ``P(|C(\tau_j)|, \sigma)`` is either a Gamma or Truncated Normal distribution with mean and standard deviation
given by ``|C(\tau_j)|`` and ``\sigma`` respectively.

The support for the Truncated Normal and Gamma distributions is ``[0,\infty)``.
Note that in the case that a Normal distribution is used it is possible for ``C_{\rm noisy}(\tau_i)``
to have a different sign that ``C(\tau_i)``, where this is not possible with the other two distributions.
"""
function add_noise(;
    # KEYWORD ARGUMENTS
    Cτ_exact::AbstractVector{T},
    τ::AbstractVector{T},
    σ::T,
    ξ::T,
    sum_rule::Function = C0 -> 1 - C0,
    noise_type::Symbol = :TruncatedNormal
) where {T<:AbstractFloat}

    Cτ_noisy = zero(Cτ_exact)
    add_noise!(
        Cτ_noisy,
        Cτ_exact = Cτ_exact,
        τ = τ,
        σ = σ,
        ξ = ξ,
        sum_rule = sum_rule = sum_rule,
        noise_type = noise_type
    )

    return Cτ_noisy
end

function add_noise(
    # ARGUMENTS
    N_samples::Int;
    # KEYWORD ARGUMENTS
    Cτ_exact::AbstractVector{T},
    τ::AbstractVector{T},
    σ::T,
    ξ::T,
    sum_rule::Function = C0 -> 1 - C0,
    noise_type::Symbol = :TruncatedNormal
) where {T<:AbstractFloat}

    Cτ_noisy = zeros(T, length(τ), N_samples)
    for column in eachcol
        add_noise!(
            column,
            Cτ_exact = Cτ_exact,
            τ = τ,
            σ = σ,
            ξ = ξ,
            sum_rule = sum_rule,
            noise_type = noise_type
        )
    end

    return Cτ_noisy
end


@doc raw"""
    add_noise!(
        # ARGUMENTS
        Cτ_noisy::AbstractVector{T};
        # KEYWORD ARGUMENTS
        Cτ_exact::AbstractVector{T},
        τ::AbstractVector{T},
        σ::T,
        ξ::T,
        sum_rule::Function = C0 -> 1 - C0,
        noise_type::Symbol = :TruncatedNormal
    ) where {T<:AbstractFloat}

    add_noise!(
        # ARGUMENTS
        Cτ_noisy::AbstractMatrix{T};
        # KEYWORD ARGUMENTS
        Cτ_exact::AbstractVector{T},
        τ::AbstractVector{T},
        σ::T,
        ξ::T,
        sum_rule::Function = C0 -> 1 - C0,
        noise_type::Symbol = :TruncatedNormal
    ) where {T<:AbstractFloat}

Add noise to an imaginary time correlation function ``C(\tau)`` that is exponentially correlated in imaginary time.

## Arguments

- `Cτ_noisy`: The array to which the noisy ``C(\tau)`` data will be written. If a matrix and not a vector then the columns correspond to samples and the rows to the imaginary time slices ``\tau``.

## Keyword Arguments

- `Cτ_exact::AbstractVector{T}`: Vector containing the exact values for ``C(\tau)``.
- `τ::AbstractVector{T}`: Vector specifying the imaginary time ``\tau`` grid that ``C(\tau)`` is evaluated on. Assumes that the last element equals the inverse temperature, i.e. `τ[end] = β`.
- `σ::T`: Standard deviation of the noise; controls the typical amplitude of the error.
- `ξ::T`: Correlation length associated with the noise in imaginary time.
- `sum_rule::Function = C0 -> 1 - C0`: Enforces sum rule, or bounday condition in imaginay time. Default behavior assumes a fermionic correlation function, enforcing that $C(\beta) = 1 - C(0)$.
- `noise_type::Symbol = :TruncatedNormal`: Distribution that the noise is sampled from prior to correlations being introduced in imaginary time. Available options are `noise_type ∈ (:TruncatedNormal, :Gamma, :Normal)`.
"""
function add_noise!(
    # ARGUMENTS
    Cτ_noisy::AbstractVector{T};
    # KEYWORD ARGUMENTS
    Cτ_exact::AbstractVector{T},
    τ::AbstractVector{T},
    σ::T,
    ξ::T,
    sum_rule::Function = C0 -> 1 - C0,
    noise_type::Symbol = :TruncatedNormal
) where {T<:AbstractFloat}

    @assert noise_type ∈ (:TruncatedNormal, :Gamma, :Normal)
    @assert size(Cτ_noisy, 1) == length(τ) == length(Cτ_exact)
    
    if noise_type == :TruncatedNormal
        R = generate_truncated_normal_noise(Cτ_exact, σ)
        @. R = sign(Cτ_exact) * R - Cτ_exact
    elseif noise_type == :Gamma
        R = generate_gamma_noise(Cτ_exact, σ)
        @. R = sign(Cτ_exact) * R - Cτ_exact
    elseif noise_type == :Normal
        R = σ * randn(length(Cτ_exact))
    end

    _add_noise!(Cτ_noisy, Cτ_exact, τ, R, ξ, sum_rule)

    return nothing
end

function add_noise!(
    # ARGUMENTS
    Cτ_noisy::AbstractMatrix{T};
    # KEYWORD ARGUMENTS
    Cτ_exact::AbstractVector{T},
    τ::AbstractVector{T},
    σ::T,
    ξ::T,
    sum_rule::Function = C0 -> 1 - C0,
    noise_type::Symbol = :TruncatedNormal
) where {T<:AbstractFloat}

    @assert noise_type ∈ (:TruncatedNormal, :Gamma, :Normal)
    @assert length(τ) == length(Cτ_exact) == length(Cτ_noisy)

    for column in eachcol(Cτ_noisy)
        add_noise!(
            column,
            Cτ_exact = Cτ_exact,
            τ = τ,
            σ = σ,
            ξ = ξ,
            sum_rule = sum_rule,
            noise_type = noise_type
        )
    end

    return nothing
end


# calculate Cτ_noisy = Cτ_exact + R∘exp(-τ/ξ) where ∘ denotes a convolution
function _add_noise!(
        Cτ_noisy::AbstractVector{T},
        Cτ_exact::AbstractVector{T},
        τ::AbstractVector{T},
        R::AbstractVector{T},
        ξ::T,
        sum_rule::Function
) where {T<:AbstractFloat}

    # length of imaginary-time axis
    Lτ = length(τ) - 1

    # get the inverse temperature
    β = τ[end]

    # R and τ arrays on interval τ ∈[0, β-Δτ]
    τ′ = @view τ[1:Lτ]
    R′ = @view R[1:Lτ]

    # calculate exponential smoothing function
    f  = @. exp(-min(τ′, β-τ′)/ξ)
    V  = sqrt(sum(fi^2 for fi in f))
    f /= V

    # convolve iid noise with exponential to generate correlated noise
    Cτ_noisy[1:Lτ] .= real.(ifft( fft(R′) .* fft(f) ))

    # generate noisy correlation data by summing true correlation with the correlated noise
    @. Cτ_noisy = Cτ_exact + Cτ_noisy

    # apply sum rule for τ = β
    Cτ_noisy[end] = sum_rule(Cτ_noisy[1])

    return nothing
end

# generate gamma distributed random noise
function generate_gamma_noise(μ::AbstractVector{T}, σ::T) where {T<:AbstractFloat}

    R = zero(μ)
    for i in eachindex(R)
        α = (abs(μ[i])/σ)^2
        θ = σ^2/abs(μ[i])
        dist = Gamma(α, θ)
        R[i] = rand(dist)
    end

    return R
end