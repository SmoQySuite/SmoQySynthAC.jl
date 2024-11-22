var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/#Spectral-to-Correlation-Function","page":"API","title":"Spectral to Correlation Function","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"spectral_to_imaginary_time_correlation_function\nspectral_to_matsubara_correlation_function","category":"page"},{"location":"api/","page":"API","title":"API","text":"spectral_to_imaginary_time_correlation_function\nspectral_to_matsubara_correlation_function","category":"page"},{"location":"api/#SmoQySynthAC.spectral_to_imaginary_time_correlation_function","page":"API","title":"SmoQySynthAC.spectral_to_imaginary_time_correlation_function","text":"spectral_to_imaginary_time_correlation_function(;\n    # KEYWORD ARGUMENTS\n    τ::AbstractVector{T},\n    β::T,\n    spectral_function::Function,\n    kernel_function::Function,\n    tol::T = 1e-10,\n) where {T<:AbstractFloat}\n\nCalculate and return the imaginary-time correlation function\n\nC(tau) = int_-infty^infty domega  K_beta(omega tau)  A(omega)\n\non a grid of tau (τ) values, given a spectral function A(omega) (spectral_function) and kernel function K_beta(omegatau) (kernel_function). This integral is evaluated within a specified tolerance tol.\n\nArguments\n\nτ::AbstractVector{T}: Vector of imaginary time such that τ[end] = β equal the inverse temperature.\nspectral_function::Function: The spectral function A(omega) that takes a single argument.\nkernel_function::Function: The kernel function K_beta(omegatau) that takes three arguments as shown.\ntol::T = 1e-10: Specified precision with which C(tau) is evaluated.\n\n\n\n\n\n","category":"function"},{"location":"api/#SmoQySynthAC.spectral_to_matsubara_correlation_function","page":"API","title":"SmoQySynthAC.spectral_to_matsubara_correlation_function","text":"spectral_to_matsubara_correlation_function(;\n    # KEYWORD ARGUMENTS\n    n::AbstractVector{Int},\n    β::T,\n    spectral_function::Function,\n    kernel_function::Function,\n    tol::T = 1e-10,\n) where {T<:AbstractFloat}\n\nCalculate and return the Matsubara correlation function\n\nC(rm i omega_n) = int_-infty^infty domega  K_beta(omega omega_n)  A(omega)\n\nfor a vector of n in mathbbZ values, given a spectral function A(omega) (spectral_function) and kernel function K_beta(omegaomega_n) (kernel_function). This integral is evaluated within a specified tolerance tol.\n\nNote that the kernel function should be called as kernel_function(ω, n, β) where n specifies the Matsubara frequency, which is evaluated internally as either omega_n = (2n+1)pibeta or omega_n = 2npibeta depending on whether the kernel function is fermionic of bosonic respectively.\n\nArguments\n\nn::AbstractVector{Int}: Vector of integers specifying Matsubara frequencies for which C(rm iomega_n) will be evaluated.\nspectral_function::Function: The spectral function A(omega) that takes a single argument.\nkernel_function::Function: The kernel function K_beta(omegaomega_n) that takes three arguments as shown.\ntol::T = 1e-10: Specified precision with which C(rm iomega_n) is evaluated.\n\n\n\n\n\n","category":"function"},{"location":"api/#Add-Noise-to-an-Imaginary-Time-Correlation-Function","page":"API","title":"Add Noise to an Imaginary Time Correlation Function","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"add_noise\nadd_noise!","category":"page"},{"location":"api/","page":"API","title":"API","text":"add_noise\nadd_noise!","category":"page"},{"location":"api/#SmoQySynthAC.add_noise","page":"API","title":"SmoQySynthAC.add_noise","text":"add_noise(;\n    # KEYWORD ARGUMENTS\n    Cτ_exact::AbstractVector{T},\n    τ::AbstractVector{T},\n    σ::T,\n    ξ::T,\n    sum_rule::Function = C0 -> 1 - C0,\n    noise_type::Symbol = :TruncatedNormal\n) where {T<:AbstractFloat}\n\nadd_noise(\n    # ARGUMENTS\n    N_samples::Int;\n    # KEYWORD ARGUMENTS\n    Cτ_exact::AbstractVector{T},\n    τ::AbstractVector{T},\n    σ::T,\n    ξ::T,\n    sum_rule::Function = C0 -> 1 - C0,\n    noise_type::Symbol = :TruncatedNormal\n) where {T<:AbstractFloat}\n\nAdd noise to an imaginary time correlation function C(tau) that is exponentially correlated in imaginary time.\n\nArguments\n\nN_samples (optional): Number of random samples to generate. If passed this function returns a matrix where the columsn correspond to the different samples.\n\nKeyword Arguments\n\nCτ_exact::AbstractVector{T}: Vector containing the exact values for C(tau).\nτ::AbstractVector{T}: Vector specifying the imaginary time tau grid that C(tau) is evaluated on. Assumes that the last element equals the inverse temperature, i.e. τ[end] = β.\nσ::T: Standard deviation of the noise; controls the typical amplitude of the error.\nξ::T: Correlation length associated with the noise in imaginary time.\nsum_rule::Function = C0 -> 1 - C0: Enforces sum rule, or bounday condition in imaginay time. Default behavior assumes a fermionic correlation function, enforcing that C(beta) = 1 - C(0).\nnoise_type::Symbol = :TruncatedNormal: Distribution that the noise is sampled from prior to correlations being introduced in imaginary time. Available options are noise_type ∈ (:TruncatedNormal, :Gamma, :Normal).\n\nAdditional Information\n\nIf the correlation function C(tau_i) the corresponding noisy correlation function is given by\n\nC_rm noisy(tau_i) = C(tau_i) + fracsum_j e^-tau_j-tau_ixi R_jsum_j e^-2tau_j-tau_ixi\n\nwhere the sum is performed assuming periodic boundary conditions. If noise_type = :Normal then the random numbers of sampled according to R_j sim rm Normal(0sigma). Otherwise,\n\nR_j sim rm sign(C(tau_j)) cdot P(C(tau_j) sigma) - C(tau_j)\n\nwhere P(C(tau_j) sigma) is either a Gamma or Truncated Normal distribution with mean and standard deviation given by C(tau_j) and sigma respectively.\n\nThe support for the Truncated Normal and Gamma distributions is 0infty). Note that in the case that a Normal distribution is used it is possible for C_rm noisy(tau_i) to have a different sign that C(tau_i), where this is not possible with the other two distributions.\n\n\n\n\n\n","category":"function"},{"location":"api/#SmoQySynthAC.add_noise!","page":"API","title":"SmoQySynthAC.add_noise!","text":"add_noise!(\n    # ARGUMENTS\n    Cτ_noisy::AbstractVector{T};\n    # KEYWORD ARGUMENTS\n    Cτ_exact::AbstractVector{T},\n    τ::AbstractVector{T},\n    σ::T,\n    ξ::T,\n    sum_rule::Function = C0 -> 1 - C0,\n    noise_type::Symbol = :TruncatedNormal\n) where {T<:AbstractFloat}\n\nadd_noise!(\n    # ARGUMENTS\n    Cτ_noisy::AbstractMatrix{T};\n    # KEYWORD ARGUMENTS\n    Cτ_exact::AbstractVector{T},\n    τ::AbstractVector{T},\n    σ::T,\n    ξ::T,\n    sum_rule::Function = C0 -> 1 - C0,\n    noise_type::Symbol = :TruncatedNormal\n) where {T<:AbstractFloat}\n\nAdd noise to an imaginary time correlation function C(tau) that is exponentially correlated in imaginary time.\n\nArguments\n\nCτ_noisy: The array to which the noisy C(tau) data will be written. If a matrix and not a vector then the columns correspond to samples and the rows to the imaginary time slices tau.\n\nKeyword Arguments\n\nCτ_exact::AbstractVector{T}: Vector containing the exact values for C(tau).\nτ::AbstractVector{T}: Vector specifying the imaginary time tau grid that C(tau) is evaluated on. Assumes that the last element equals the inverse temperature, i.e. τ[end] = β.\nσ::T: Standard deviation of the noise; controls the typical amplitude of the error.\nξ::T: Correlation length associated with the noise in imaginary time.\nsum_rule::Function = C0 -> 1 - C0: Enforces sum rule, or bounday condition in imaginay time. Default behavior assumes a fermionic correlation function, enforcing that C(beta) = 1 - C(0).\nnoise_type::Symbol = :TruncatedNormal: Distribution that the noise is sampled from prior to correlations being introduced in imaginary time. Available options are noise_type ∈ (:TruncatedNormal, :Gamma, :Normal).\n\n\n\n\n\n","category":"function"},{"location":"api/#Kernel-Functions","page":"API","title":"Kernel Functions","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"kernel_tau_fermi\nkernel_tau_bose\nkernel_tau_sym_bose\nkernel_mat_fermi\nkernel_mat_bose","category":"page"},{"location":"api/","page":"API","title":"API","text":"kernel_tau_fermi\nkernel_tau_bose\nkernel_tau_sym_bose\nkernel_mat_fermi\nkernel_mat_bose","category":"page"},{"location":"api/#SmoQySynthAC.kernel_tau_fermi","page":"API","title":"SmoQySynthAC.kernel_tau_fermi","text":"kernel_tau_fermi(ω::T, τ::T, β::T) where {T<:AbstractFloat}\n\nThe imaginary time fermionic kernel\n\nbeginalign\nK_beta(omegatau)  = overbraceleft(frace^-tauomega1+e^-betaomegaright)^textnumerically unstable \n                      = underbraceleft( e^tauomega + e^(tau-beta)omega right)^-1_textnumerically stable\nendalign\n\nwhere it is assumed that tau in 0beta).\n\n\n\n\n\n","category":"function"},{"location":"api/#SmoQySynthAC.kernel_tau_bose","page":"API","title":"SmoQySynthAC.kernel_tau_bose","text":"kernel_tau_bose(ω::T, τ::T, β::T) where {T<:AbstractFloat}\n\nThe imaginary time bosonic kernel\n\nbeginalign\nK_beta(omegatau)  = overbraceleft(frace^-tauomega1-e^-betaomegaright)^textnumerically unstable \n                      = underbraceleft( e^tauomega - e^(tau-beta)omega right)^-1_textnumerically stable\nendalign\n\nwhere it is assumed that tau in 0beta).\n\n\n\n\n\n","category":"function"},{"location":"api/#SmoQySynthAC.kernel_tau_sym_bose","page":"API","title":"SmoQySynthAC.kernel_tau_sym_bose","text":"kernel_tau_sym_bose(ω::T, τ::T, β::T) where {T<:AbstractFloat}\n\nThe imaginary time symmetrized bosonic kernel\n\nK_beta(omegatau) = frace^-tauomega + e^-(beta-tau)omega1-e^-betaomega\n\nwhere it is assumed that tau in 0beta).\n\n\n\n\n\n","category":"function"},{"location":"api/#SmoQySynthAC.kernel_mat_fermi","page":"API","title":"SmoQySynthAC.kernel_mat_fermi","text":"kernel_mat_fermi(ω::T, n::Int, β::T) where {T<:AbstractFloat}\n\nThe fermionic matsubara frequency kernel\n\nK_beta(omega omega_n) = frac1rm iomega_n - omega\n\nwhere omega_n = (2n+1)pibeta for fermions with n in mathbbZ.\n\n\n\n\n\n","category":"function"},{"location":"api/#SmoQySynthAC.kernel_mat_bose","page":"API","title":"SmoQySynthAC.kernel_mat_bose","text":"kernel_mat_bose(ω::T, n::Int, β::T) where {T<:AbstractFloat}\n\nThe bosonic matsubara frequency kernel\n\nK_beta(omega omega_n) = frac1rm iomega_n - omega\n\nwhere omega_n = 2npibeta for bosons with n in mathbbZ.\n\n\n\n\n\n","category":"function"},{"location":"api/#Utilitie-Functions","page":"API","title":"Utilitie Functions","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"fermi\nbose","category":"page"},{"location":"api/","page":"API","title":"API","text":"fermi\nbose","category":"page"},{"location":"api/#SmoQySynthAC.fermi","page":"API","title":"SmoQySynthAC.fermi","text":"fermi(ϵ::T, β::T) where {T<:AbstractFloat}\n\nThe Fermi-Dirac function\n\nbeginalign\nf_beta(epsilon)  = overbraceleft( frac1e^betaepsilon + 1right)^textnumerically unstable \n                   = underbracefrac12left( 1 - tanhleft(fracbetaepsilon2right) right)_textnumerically stable\nendalign\n\nwhere epsilon is energy and beta is inverse temperature.\n\n\n\n\n\n","category":"function"},{"location":"api/#SmoQySynthAC.bose","page":"API","title":"SmoQySynthAC.bose","text":"bose(ϵ::T, β::T) where {T<:AbstractFloat}\n\nThe Bose-Einstein function\n\nbeginalign\nn_beta(epsilon)  = overbraceleft( frac1e^betaepsilon - 1right)^textnumerically unstable \n                   = underbracefrac12left(cothleft(fracbetaepsilon2right) - 1 right)_textnumerically stable\nendalign\n\nwhere epsilon is energy and beta is inverse temperature.\n\n\n\n\n\n","category":"function"},{"location":"usage/","page":"Usage","title":"Usage","text":"EditURL = \"../usage.jl\"","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Here we give a basic example demonstrating some of the functionality of the SmoQySynthAC.jl package.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"using SmoQySynthAC\nusing Distributions\nusing CairoMakie\n\n# render figures using SVG backend so they are sharp\nCairoMakie.activate!(type = \"svg\")","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"In this example we will work with the single-particle imaginary time fermion Green's function which is given by","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"G(tau) = int_-infty^infty K_beta(omegatau) A(omega)","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"where A(omega) is the spectral function and","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"K_beta(omegatau) = frace^-tau omega1 + e^-beta omega","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"is the kernel function where beta = 1T is the inverse temperature and it is assumed that tau in 0 beta).","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"As a first step in demonstrating the functionality of SmoQySynthAC.jl package, let us define a synthetic spectral function A(omega). For convenience we will do this using the Distributions.jl package. Here we define the spectral function as Lorentzian (Cauchy) distribution in between two Normal distributions.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"# define spectral function distribution\nspectral_dist = MixtureModel(\n    [Normal(-2.0,0.7), Cauchy(0.0, 0.3), Normal(+2.0,0.7)],\n    [0.2, 0.4, 0.4]\n)\n\n# define method to evaluate spectral function\nspectral_function = ω -> pdf(spectral_dist, ω)","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Now let us quickly plot our spectral function so we can see what it looks like.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"ωmin = -6.0\nωmax = 6.0\nω = collect(range(start = ωmin, stop = ωmax, length = 1000))\nA = spectral_function.(ω)\n\nfig = Figure(\n    size = (700, 400),\n    fonts = (; regular= \"CMU Serif\"),\n    figure_padding = 10\n)\n\nax = Axis(\n    fig[1, 1],\n    aspect = 7/4,\n    xlabel = L\"\\omega\", ylabel = L\"A(\\omega)\",\n    xlabelsize = 30, ylabelsize = 30,\n    xticklabelsize = 24, yticklabelsize = 24,\n)\n\nlines!(\n    ω, A,\n    linewidth = 2,\n    alpha = 1.0,\n    color = :black,\n    linestyle = :solid\n)\n\nxlims!(ax, ωmin, ωmax)\nylims!(ax, 0.0, 1.05*maximum(A))\n\nfig","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"The next step is to define the inverse temperature beta, discretization in imaginary Deltatau and corresponding tau grid.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"# Set inverse temperature.\nβ = 10.0\n\n# Set discretization in imaginary time.\nΔτ = 0.05\n\n# Calculate corresponding imaginary time grid.\nτ = collect(range(start = 0.0, stop = β, step = Δτ));\nnothing #hide","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Now we can calculate G(tau) using the spectral_to_imaginary_time_correlation_function method and appropriate kernel function kernel_tau_fermi.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"# Calculate imaginary time Green's function.\nGτ = spectral_to_imaginary_time_correlation_function(\n    τ = τ,\n    β = β,\n    spectral_function = spectral_function,\n    kernel_function = kernel_tau_fermi,\n    tol = 1e-10\n);\nnothing #hide","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"We can similary calculate the Matsubara Green's function G(textiomega_n) using the function spectral_to_matsubara_correlation_function function with the kernel function kernel_mat_fermi.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"# Define Matsubara frequency grid in terms of integers n where ωₙ = (2n+1)π/β.\nn = collect(-250:250)\n\n# Calculate Matsubara Green's function.\nGn = spectral_to_matsubara_correlation_function(;\n    n = n,\n    β = β,\n    spectral_function = spectral_function,\n    kernel_function = kernel_mat_fermi,\n    tol = 1e-10,\n);\nnothing #hide","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"The resulting real and imaginary parts of G(textiomega_n) are plotted below.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"ωn = @. (2*n+1)*π/β\n\nfig = Figure(\n    size = (700, 500),\n    fonts = (; regular= \"CMU Serif\"),\n    figure_padding = 10\n)\n\nax = Axis(fig[1, 1],\n    aspect = 7/5,\n    xlabel = L\"\\omega_n\",\n    ylabel = L\"G_\\sigma(\\text{i}\\omega_n)\",\n    xlabelsize = 30, ylabelsize = 30,\n    xticklabelsize = 24, yticklabelsize = 24,\n)\n\nxlims!(ax, minimum(ωn), maximum(ωn))\n\nlines!(\n    ωn, real.(Gn),\n    linewidth = 2, alpha = 2.0, color = :red, linestyle = :solid, label = L\"\\text{Re}[G_\\sigma(\\text{i}\\omega_n)]\"\n)\n\nlines!(\n    ωn, imag.(Gn),\n    linewidth = 2, alpha = 1.5, color = :green, linestyle = :solid, label = L\"\\text{Im}[G_\\sigma(\\text{i}\\omega_n)]\"\n)\n\naxislegend(\n    ax, halign = :left, valign = :top, labelsize = 30\n)\n\nfig","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Having calculated the exact G(tau) function, let us now add some noise to it using the add_noise method.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Gτ_noisy = add_noise(\n    Cτ_exact = Gτ,\n    τ = τ,\n    σ = 1.0e-2,\n    ξ = 0.5,\n    sum_rule = C0 -> 1 - C0, # enforces fermionic boundary condition that C(τ=β) = 1 - C(τ=0)\n    noise_type = :TruncatedNormal\n);\nnothing #hide","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Let us now plot what both G(tau) and G_rm noisy(tau).","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"fig = Figure(\n    size = (700, 500),\n    fonts = (; regular= \"CMU Serif\"),\n    figure_padding = 10\n)\n\nax = Axis(fig[1, 1],\n    aspect = 7/5,\n    xlabel = L\"\\tau\",\n    ylabel = L\"G_\\sigma(\\tau)\",\n    xlabelsize = 30, ylabelsize = 30,\n    xticklabelsize = 24, yticklabelsize = 24,\n)\n\nxlims!(ax, 0.0, β)\nylims!(ax, 0.0, 1.0)\n\nlines!(\n    τ, Gτ,\n    linewidth = 2, alpha = 2.0, color = :black, linestyle = :solid, label = \"Exact\"\n)\n\nlines!(\n    τ, Gτ_noisy,\n    linewidth = 2, alpha = 1.5, color = :red, linestyle = :solid, label = \"Noisy\"\n)\n\naxislegend(\n    ax, halign = :left, valign = :top, labelsize = 30\n)\n\nfig","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = SmoQySynthAC","category":"page"},{"location":"#SmoQySynthAC","page":"Home","title":"SmoQySynthAC","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"SmoQySynthAC.jl is a package developed to generate synthetic imaginary time correlation data mean to resemble noisy quantum Monte Carlo (QMC) measurements, where the introduced noise is correlated in imaginary time. This synthetic data can then be used to benchmark analytic continuation (AC) methods, or numerical methods for transforming noisy imaginary time correlation data to Matsubara space. The approach used to generate the noisy synthetic data is based method introduced in Shao et. al., the citation for which is given below.","category":"page"},{"location":"","page":"Home","title":"Home","text":"@article{PhysRevX.7.041072,\n  title = {Nearly Deconfined Spinon Excitations in the Square-Lattice Spin-$1/2$ Heisenberg Antiferromagnet},\n  author = {Shao, Hui and Qin, Yan Qi and Capponi, Sylvain and Chesi, Stefano and Meng, Zi Yang and Sandvik, Anders W.},\n  journal = {Phys. Rev. X},\n  volume = {7},\n  issue = {4},\n  pages = {041072},\n  numpages = {26},\n  year = {2017},\n  month = {Dec},\n  publisher = {American Physical Society},\n  doi = {10.1103/PhysRevX.7.041072},\n  url = {https://link.aps.org/doi/10.1103/PhysRevX.7.041072}\n}","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install SmoQySynthAC.jl, simply open the Julia REPL and run the commands","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> ]\npkg> add SmoQySynthAC","category":"page"},{"location":"","page":"Home","title":"Home","text":"or equivalently via Pkg do","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg; Pkg.add(\"SmoQySynthAC\")","category":"page"},{"location":"#Funding","page":"Home","title":"Funding","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The development of this code was supported by the U.S. Department of Energy, Office of Science, Basic Energy Sciences, under Award Number DE-SC0022311.","category":"page"}]
}
