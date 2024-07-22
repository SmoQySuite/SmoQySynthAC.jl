# Here we give a basic example demonstrating some of the functionality of the
# [SmoQySynthAC.jl](https://github.com/SmoQySuite/SmoQySynthAC.jl.git) package.

using SmoQySynthAC
using Distributions
using CairoMakie

## render figures using SVG backend so they are sharp
CairoMakie.activate!(type = "svg")

# In this example we will work with the single-particle imaginary time fermion Green's function
# which is given by
# ```math
# G(\tau) = \int_{-\infty}^\infty K_\beta(\omega,\tau) A(\omega)
# ```
# where ``A(\omega)`` is the spectral function and
# ```math
# K_\beta(\omega,\tau) = \frac{e^{-\tau \omega}}{1 + e^{-\beta \omega}}
# ```
# is the kernel function where ``\beta = 1/T`` is the inverse temperature and it is assumed that
# ``\tau \in [0, \beta)``.
#
# As a first step in demonstrating the functionality of [SmoQySynthAC.jl](https://github.com/SmoQySuite/SmoQySynthAC.jl.git) package,
# let us define a synthetic spectral function ``A(\omega)``.
# For convenience we will do this using the [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl.git) package.
# We will define a spectral function with a Lorentzian (Cauchy) distribution in centered between two Normal distributions on either side.

## define spectral function distribution
spectral_dist = MixtureModel(
    [Normal(-2.0,0.7), Cauchy(0.0, 0.3), Normal(+2.0,0.7)],
    [0.2, 0.4, 0.4]
)

## define method to evaluate spectral function
spectral_function = ω -> pdf(spectral_dist, ω)

# Now let us quickly plot our spectral function so we can see what it looks like.

ωmin = -6.0
ωmax = 6.0
ω = collect(range(start = ωmin, stop = ωmax, length = 1000))
A = spectral_function.(ω)

fig = Figure(
    size = (700, 400),
    fonts = (; regular= "CMU Serif"),
    figure_padding = 10
)

ax = Axis(
    fig[1, 1],
    aspect = 7/4,
    xlabel = L"\omega", ylabel = L"A(\omega)",
    xlabelsize = 30, ylabelsize = 30,
    xticklabelsize = 24, yticklabelsize = 24,
)

lines!(
    ω, A,
    linewidth = 2,
    alpha = 1.0,
    color = :black,
    linestyle = :solid
)

xlims!(ax, ωmin, ωmax)
ylims!(ax, 0.0, 1.05*maximum(A))

fig

# The next step is to define the inverse temperature ``\beta``, discretization in imaginary ``\Delta\tau``
# and corresponding ``\tau`` grid.

## Set inverse temperature.
β = 10.0

## Set discretization in imaginary time.
Δτ = 0.05

## Calculate corresponding imaginary time grid.
τ = collect(range(start = 0.0, stop = β, step = Δτ));

# Now we can calculate ``G(\tau)`` using the [`spectral_to_imaginary_time_correlation_function`](@ref)
# method and appropriate kernel funciton [`kernel_tau_fermi`](@ref).

## Calculate imaginary time Green's function.
Gτ = spectral_to_imaginary_time_correlation_function(
    τ = τ,
    β = β,
    spectral_function = spectral_function,
    kernel_function = kernel_tau_fermi,
    tol = 1e-10
);

# We can similary calculate the Matsubara Green's function ``G(\text{i}\omega_n)``
# using the function [`spectral_to_matsubara_correlation_function`](@ref) function
# with the kernel function [`kernel_mat_fermi`](@ref).

## Define Matsubara frequency grid in terms of integers n where ωₙ = (2n+1)π/β.
n = collect(-250:250)

## Calculate Matsubara Green's function.
Gn = spectral_to_matsubara_correlation_function(;
    n = n,
    β = β,
    spectral_function = spectral_function,
    kernel_function = kernel_mat_fermi,
    tol= 1e-10,
);

# The resulting real and imaginary parts of ``G(\text{i}\omega_n)`` are plotted below.

ωn = @. (2*n+1)*π/β

fig = Figure(
    size = (700, 500),
    fonts = (; regular= "CMU Serif"),
    figure_padding = 10
)

ax = Axis(fig[1, 1],
    aspect = 7/5,
    xlabel = L"\omega_n",
    ylabel = L"G_\sigma(\text{i}\omega_n)",
    xlabelsize = 30, ylabelsize = 30,
    xticklabelsize = 24, yticklabelsize = 24,
)

xlims!(ax, minimum(ωn), maximum(ωn))

lines!(
    ωn, real.(Gn),
    linewidth = 2, alpha = 2.0, color = :red, linestyle = :solid, label = L"\text{Re}[G_\sigma(\text{i}\omega_n)]"
)

lines!(
    ωn, imag.(Gn),
    linewidth = 2, alpha = 1.5, color = :green, linestyle = :solid, label = L"\text{Im}[G_\sigma(\text{i}\omega_n)]"
)

axislegend(
    ax, halign = :left, valign = :top, labelsize = 30
)

fig

# Having calculated the exact ``G(\tau)`` function, let us now add some noise to it using the
# [`add_noise`](@ref) method.

Gτ_noisy = add_noise(
    Cτ_exact = Gτ,
    τ = τ,
    σ = 1.0e-2,
    ξ = 0.5,
    sum_rule = C0 -> 1 - C0, # enforces fermionic boundary condition that C(τ=β) = 1 - C(τ=0)
    noise_type = :TruncatedNormal
);

# Let us now plot what both ``G(\tau)`` and ``G_{\rm noisy}(\tau)``.

fig = Figure(
    size = (700, 500),
    fonts = (; regular= "CMU Serif"),
    figure_padding = 10
)

ax = Axis(fig[1, 1],
    aspect = 7/5,
    xlabel = L"\tau",
    ylabel = L"G_\sigma(\tau)",
    xlabelsize = 30, ylabelsize = 30,
    xticklabelsize = 24, yticklabelsize = 24,
)

xlims!(ax, 0.0, β)
ylims!(ax, 0.0, 1.0)

lines!(
    τ, Gτ,
    linewidth = 2, alpha = 2.0, color = :black, linestyle = :solid, label = "Exact"
)

lines!(
    τ, Gτ_noisy,
    linewidth = 2, alpha = 1.5, color = :red, linestyle = :solid, label = "Noisy"
)

axislegend(
    ax, halign = :left, valign = :top, labelsize = 30
)

fig