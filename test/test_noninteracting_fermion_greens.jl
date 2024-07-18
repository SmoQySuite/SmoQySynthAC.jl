@testitem "Test Non-Interacting Imaginary-Time Green's Function" begin

    # function to define spectral functions in terms of simple pole expansion
    simple_pole(z; location, residue) = residue/(z-location)
    simple_poles(z; locations, residues) = sum(simple_pole(z, location = locations[i], residue = residues[i]) for i in eachindex(locations))
    spe_spectral_function(ω; locations, residues, η=1.0e-9) = -inv(π)*imag(simple_poles(ω + η*im, locations = locations, residues = residues))

    # non-interacting single-particle fermion green's function
    noninteracting_greens(τ, β, ϵ) = inv(exp(τ*ϵ) + exp((τ-β)*ϵ))

    # single-particle state energy
    ϵ = 0.2

    # inverse temperature
    β = 3.0

    # discretization in imaginary time
    Δτ = 0.10

    # imaginary time grid
    τ = collect(range(start = 0.0, stop = β, step = Δτ))

    # exact non-interacting greens
    G_exact = noninteracting_greens.(τ, β, ϵ)

    # numerical fudge factor for pole-represention of non-interacting greens
    η = 1.0e-8

    # define spectral function
    locations = [ϵ - η*im]
    residues = [1.0 + 0.0*im]
    spectral_function = ω -> spe_spectral_function(ω, locations = locations, residues = residues, η = η)

    # calculate non-interacting greens by integrating the spectral function and kernel
    G_int = spectral_to_imaginary_time_correlation_function(
        τ = τ,
        β = β,
        spectral_function = spectral_function,
        kernel_function = kernel_tau_fermi,
        tol = 1e-10
    )

    # test that the two greens agree
    @test maximum(abs.(G_int - G_exact)) < sqrt(η)

end