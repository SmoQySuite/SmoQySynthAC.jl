```@meta
EditURL = "../usage.jl"
```

Here we give a basic example demonstrating some of the functionality of the
[SmoQySynthAC.jl](https://github.com/SmoQySuite/SmoQySynthAC.jl.git) package.

````@example usage
using SmoQySynthAC
using Distributions
using CairoMakie

# render figures using SVG backend so they are sharp
CairoMakie.activate!(type = "svg")
````

In this example we will work with the single-particle imaginary time fermion Green's function
which is given by
```math
G(\tau) = \int_{-\infty}^\infty K(\omega,\tau,\beta) A(\omega)
```
where ``A(\omega)`` is the spectral function and
```math
K(\omega,\tau,\beta) = \frac{e^{-\tau \omega}}{1 + e^{-\beta \omega}}
```
is the kernel function where ``\beta = 1/T`` is the inverse temperature and it is assumed that
``\tau \in [0, \beta)``.

As a first step in demonstrating the functionality of [SmoQySynthAC.jl](https://github.com/SmoQySuite/SmoQySynthAC.jl.git) package,
let us define a synthetic spectral function ``A(\omega)``.
For convenience we will do this using the [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl.git) package.
We will define a spectral function with a Lorentzian (Cauchy) distribution in centered between two Normal distributions on either side.

````@example usage
# define spectral distribution
spectral_dist = MixtureModel(
    [Normal(-2.0,0.7), Cauchy(0.0, 0.3), Normal(+2.0,0.7)],
    [0.2, 0.4, 0.4]
)

# define function to evaluate the spectral funciton
spectral_function = ω -> pdf(spectral_dist, ω)
````

Now let us quickly plot our spectral function so we can see what it looks like.

````@example usage
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
````

The next step is to define the inverse temperature ``\beta``, discretization in imaginary ``\Delta\tau``
and corresponding ``\tau`` grid.

````@example usage
β = 10.0
Δτ = 0.05
τ = collect(range(start = 0.0, stop = β, step = Δτ));
nothing #hide
````

Now we can calculate ``G(\tau)`` using the [`spectral_to_imaginary_time_correlation_function`](@ref)
method and appropriate kernel funciton [`kernel_tau_fermi`](@ref).

````@example usage
Gτ = spectral_to_imaginary_time_correlation_function(
    τ = τ,
    β = β,
    spectral_function = spectral_function,
    kernel_function = kernel_tau_fermi,
    tol = 1e-10
);
nothing #hide
````

Having calculated the exact ``G(\tau)`` function, let us now add some noise to it using the
[`add_noise`](@ref) method.

````@example usage
Gτ_noisy = add_noise(
    Cτ_exact = Gτ,
    τ = τ,
    σ = 1.0e-2,
    ξ = 0.5,
    sum_rule = C0 -> 1 - C0, # enforces fermionic boundary condition that C(τ=β) = 1 - C(τ=0)
    noise_type = :TruncatedNormal
);
nothing #hide
````

Let us now plot what both ``G(\tau)`` and ``G_{\rm noisy}(\tau)``.

````@example usage
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
````

