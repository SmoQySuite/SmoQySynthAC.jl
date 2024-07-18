
# Return a vector of random numbers R, where each element is sampled from a truncated normal
# distribution defined on the interval [0,∞) with mean given by correspond element of the vector |μt|,
# all sharing the same standard deviation σt
function generate_truncated_normal_noise(μt::AbstractVector{T}, σt::T) where {T<:AbstractFloat}

    R = zero(μt)
    for i in eachindex(R)
        dist = truncated_normal(abs(μt[i]), σt)
        R[i] = rand(dist)
    end

    return R
end

# Return truncated normal distribution defined with mean and standard deviation
# (μt, σt) > 0 defined on the interval [0,∞)
function truncated_normal(μt::T, σt::T) where {T<:AbstractFloat}

    lower = [zero(T), zero(T)]
    upper = [T(Inf), T(Inf)]
    p0 = [μt, σt]
    f = params -> loss(params[1], params[2], μt, σt)
    res = optimize(f, lower, upper, p0, Fminbox(LBFGS()))
    sol = Optim.minimizer(res)
    μ = sol[1]
    σ = sol[2]

    return truncated(Normal(μ,σ), 0.0, Inf)
end

# loss function used in optimization
function loss(μ, σ, μt, σt)

    μt′ = eval_μt(μ, σ)
    σt′ = eval_σt(μ, σ)
    return (μt - μt′)^2 + (σt - σt′)^2
end

function eval_σt(μ, σ)

    α = -μ/σ
    φ = eval_φ(μ, σ)
    Z = eval_Z(μ, σ)
    φoZ = φ/Z
    return σ * sqrt( 1 + α * φoZ + φoZ^2 )
end

function eval_μt(μ, σ)

    φ = eval_φ(μ, σ)
    Z = eval_Z(μ, σ)
    return μ + σ * φ/Z
end

function eval_Z(μ, σ)
    
    return (1 + erf(μ/(sqrt(2)*σ)))/2
end

function eval_φ(μ, σ)
    
    return exp(-(μ/σ)^2/2)/sqrt(2*π)
end