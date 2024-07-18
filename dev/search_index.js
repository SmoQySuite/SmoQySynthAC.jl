var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/#Spectral-to-Correlation-Function","page":"API","title":"Spectral to Correlation Function","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"spectral_to_imaginary_time_correlation_function\nspectral_to_matsubara_correlation_function","category":"page"},{"location":"api/","page":"API","title":"API","text":"spectral_to_imaginary_time_correlation_function\nspectral_to_matsubara_correlation_function","category":"page"},{"location":"api/#SmoQySynthAC.spectral_to_imaginary_time_correlation_function","page":"API","title":"SmoQySynthAC.spectral_to_imaginary_time_correlation_function","text":"spectral_to_imaginary_time_correlation_function(;\n    # KEYWORD ARGUMENTS\n    τ::AbstractVector{T},\n    β::T,\n    spectral_function::Function,\n    kernel_function::Function,\n    tol::T = 1e-10,\n) where {T<:AbstractFloat}\n\nCalculate and return the imaginary-time correlation function\n\nC(tau) = int_-infty^infty domega  K(omega tau beta)  A(omega)\n\non a grid of tau (τ) values, given a spectral function A(omega) (spectral_function) and kernel function K(omegataubeta) (kernel_function). This integral is evaluated within a specified tolerance tol.\n\nArguments\n\nτ::AbstractVector{T}: Vector of imaginary time such that τ[end] = β equal the inverse temperature.\nspectral_function::Function: The spectral function A(omega) that takes a single argument.\nkernel_function::Function: The kernel function K(omegataubeta) that takes three arguments as shown.\ntol::T = 1e-10: Specified precision with which C(tau) is evaluated.\n\n\n\n\n\n","category":"function"},{"location":"api/#SmoQySynthAC.spectral_to_matsubara_correlation_function","page":"API","title":"SmoQySynthAC.spectral_to_matsubara_correlation_function","text":"spectral_to_matsubara_correlation_function(;\n    # KEYWORD ARGUMENTS\n    n::AbstractVector{Int},\n    β::T,\n    spectral_function::Function,\n    kernel_function::Function,\n    tol::T = 1e-10,\n) where {T<:AbstractFloat}\n\nCalculate and return the imaginary-time correlation function\n\nC(tau) = int_-infty^infty domega  K(omega omega_n beta)  A(omega)\n\non a grid of tau (τ) values, given a spectral function A(omega) (spectral_function) and kernel function K(omegataubeta) (kernel_function). This integral is evaluated within a specified tolerance tol.\n\nNote that the kernel function should be called as kernel_function(ω, n, β) where n specifies the Matsubara frequency, which is evaluated internally as either omega_n = (2n+1)pibeta or omega_n = 2npibeta depending on whether the kernel function is fermionic of bosonic respectively.\n\nArguments\n\nτ::AbstractVector{T}: Vector of imaginary time such that τ[end] = β equal the inverse temperature.\nspectral_function::Function: The spectral function A(omega) that takes a single argument.\nkernel_function::Function: The kernel function K(omegataubeta) that takes three arguments as shown.\ntol::T = 1e-10: Specified precision with which C(tau) is evaluated.\n\n\n\n\n\n","category":"function"},{"location":"api/#Add-Noise-to-an-Imaginary-Time-Correlation-Function","page":"API","title":"Add Noise to an Imaginary Time Correlation Function","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"add_noise\nadd_noise!","category":"page"},{"location":"api/","page":"API","title":"API","text":"add_noise\nadd_noise!","category":"page"},{"location":"api/#SmoQySynthAC.add_noise","page":"API","title":"SmoQySynthAC.add_noise","text":"add_noise(;\n    # KEYWORD ARGUMENTS\n    Cτ_exact::AbstractVector{T},\n    τ::AbstractVector{T},\n    σ::T,\n    ξ::T,\n    sum_rule::Function = C0 -> 1 - C0,\n    noise_type::Symbol = :TruncatedNormal\n) where {T<:AbstractFloat}\n\nadd_noise(\n    # ARGUMENTS\n    N_samples::Int;\n    # KEYWORD ARGUMENTS\n    Cτ_exact::AbstractVector{T},\n    τ::AbstractVector{T},\n    σ::T,\n    ξ::T,\n    sum_rule::Function = C0 -> 1 - C0,\n    noise_type::Symbol = :TruncatedNormal\n) where {T<:AbstractFloat}\n\nAdd noise to an imaginary time correlation function C(tau) that is exponentially correlated in imaginary time.\n\nArguments\n\nN_samples (optional): Number of random samples to generate. If passed this function returns a matrix where the columsn correspond to the different samples.\n\nKeyword Arguments\n\nCτ_exact::AbstractVector{T}: Vector containing the exact values for C(tau).\nτ::AbstractVector{T}: Vector specifying the imaginary time tau grid that C(tau) is evaluated on. Assumes that the last element equals the inverse temperature, i.e. τ[end] = β.\nσ::T: Standard deviation of the noise; controls the typical amplitude of the error.\nξ::T: Correlation length associated with the noise in imaginary time.\nsum_rule::Function = C0 -> 1 - C0: Enforces sum rule, or bounday condition in imaginay time. Default behavior assumes a fermionic correlation function, enforcing that C(beta) = 1 - C(0).\nnoise_type::Symbol = :TruncatedNormal: Distribution that the noise is sampled from prior to correlations being introduced in imaginary time. Available options are noise_type ∈ (:TruncatedNormal, :Gamma, :Normal).\n\nAdditional Information\n\nIf the correlation function C(tau_i) the corresponding noisy correlation function is given by\n\nC_rm noisy(tau_i) = C(tau_i) + fracsum_j e^-tau_j-tau_ixi R_jsum_j e^-2tau_j-tau_ixi\n\nwhere the sum is performed assuming periodic boundary conditions. If noise_type = :Normal then the random numbers of sampled according to R_j sim rm Normal(0sigma). Otherwise,\n\nR_j sim rm sign(C(tau_j)) cdot P(C(tau_j) sigma) - C(tau_j)\n\nwhere P(C(tau_j) sigma) is either a Gamma or Truncated Normal distribution with mean and standard deviation given by C(tau_j) and sigma respectively.\n\nThe support for the Truncated Normal and Gamma distributions is 0infty). Note that in the case that a Normal distribution is used it is possible for C_rm noisy(tau_i) to have a different sign that C(tau_i), where this is not possible with the other two distributions.\n\n\n\n\n\n","category":"function"},{"location":"api/#SmoQySynthAC.add_noise!","page":"API","title":"SmoQySynthAC.add_noise!","text":"add_noise!(\n    # ARGUMENTS\n    Cτ_noisy::AbstractVector{T};\n    # KEYWORD ARGUMENTS\n    Cτ_exact::AbstractVector{T},\n    τ::AbstractVector{T},\n    σ::T,\n    ξ::T,\n    sum_rule::Function = C0 -> 1 - C0,\n    noise_type::Symbol = :TruncatedNormal\n) where {T<:AbstractFloat}\n\nadd_noise!(\n    # ARGUMENTS\n    Cτ_noisy::AbstractMatrix{T};\n    # KEYWORD ARGUMENTS\n    Cτ_exact::AbstractVector{T},\n    τ::AbstractVector{T},\n    σ::T,\n    ξ::T,\n    sum_rule::Function = C0 -> 1 - C0,\n    noise_type::Symbol = :TruncatedNormal\n) where {T<:AbstractFloat}\n\nAdd noise to an imaginary time correlation function C(tau) that is exponentially correlated in imaginary time.\n\nArguments\n\nCτ_noisy: The array to which the noisy C(tau) data will be written. If a matrix and not a vector then the columns correspond to samples and the rows to the imaginary time slices tau.\n\nKeyword Arguments\n\nCτ_exact::AbstractVector{T}: Vector containing the exact values for C(tau).\nτ::AbstractVector{T}: Vector specifying the imaginary time tau grid that C(tau) is evaluated on. Assumes that the last element equals the inverse temperature, i.e. τ[end] = β.\nσ::T: Standard deviation of the noise; controls the typical amplitude of the error.\nξ::T: Correlation length associated with the noise in imaginary time.\nsum_rule::Function = C0 -> 1 - C0: Enforces sum rule, or bounday condition in imaginay time. Default behavior assumes a fermionic correlation function, enforcing that C(beta) = 1 - C(0).\nnoise_type::Symbol = :TruncatedNormal: Distribution that the noise is sampled from prior to correlations being introduced in imaginary time. Available options are noise_type ∈ (:TruncatedNormal, :Gamma, :Normal).\n\n\n\n\n\n","category":"function"},{"location":"api/#Kernel-Functions","page":"API","title":"Kernel Functions","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"kernel_tau_fermi\nkernel_tau_bose\nkernel_tau_sym_bose\nkernel_mat_fermi\nkernel_mat_bose","category":"page"},{"location":"api/","page":"API","title":"API","text":"kernel_tau_fermi\nkernel_tau_bose\nkernel_tau_sym_bose\nkernel_mat_fermi\nkernel_mat_bose","category":"page"},{"location":"api/#SmoQySynthAC.kernel_tau_fermi","page":"API","title":"SmoQySynthAC.kernel_tau_fermi","text":"kernel_tau_fermi(ω::T, τ::T, β::T) where {T<:AbstractFloat}\n\nThe imaginary time fermionic kernel\n\nK_beta(omegatau) = frace^-tauomega1+e^-betaomega\n\nwhere it is assumed that tau in 0beta).\n\n\n\n\n\n","category":"function"},{"location":"api/#SmoQySynthAC.kernel_tau_bose","page":"API","title":"SmoQySynthAC.kernel_tau_bose","text":"kernel_tau_bose(ω::T, τ::T, β::T) where {T<:AbstractFloat}\n\nThe imaginary time bosonic kernel\n\nK_beta(omegatau) = frace^-tauomega1-e^-betaomega\n\nwhere it is assumed that tau in 0beta).\n\n\n\n\n\n","category":"function"},{"location":"api/#SmoQySynthAC.kernel_tau_sym_bose","page":"API","title":"SmoQySynthAC.kernel_tau_sym_bose","text":"kernel_tau_sym_bose(ω::T, τ::T, β::T) where {T<:AbstractFloat}\n\nThe imaginary time symmetrized bosonic kernel\n\nK_beta(omegatau) = frace^-tauomega + e^-(beta-tau)omega1-e^-betaomega\n\nwhere it is assumed that tau in 0beta).\n\n\n\n\n\n","category":"function"},{"location":"api/#SmoQySynthAC.kernel_mat_fermi","page":"API","title":"SmoQySynthAC.kernel_mat_fermi","text":"kernel_mat_fermi(ω::T, n::Int, β::T) where {T<:AbstractFloat}\n\nThe fermionic matsubara frequency kernel\n\nK_beta(omega rm iomega_n) = frac1rm iomega_n - omega\n\nwhere omega_n = (2n+1)pibeta for fermions with n in mathbbZ.\n\n\n\n\n\n","category":"function"},{"location":"api/#SmoQySynthAC.kernel_mat_bose","page":"API","title":"SmoQySynthAC.kernel_mat_bose","text":"kernel_mat_bose(ω::T, n::Int, β::T) where {T<:AbstractFloat}\n\nThe bosonic matsubara frequency kernel\n\nK_beta(omega rm iomega_n) = frac1rm iomega_n - omega\n\nwhere omega_n = 2npibeta for bosons with n in mathbbZ.\n\n\n\n\n\n","category":"function"},{"location":"usage/","page":"Usage","title":"Usage","text":"EditURL = \"../usage.jl\"","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Here we give a basic example demonstrating some of the functionality of the SmoQySynthAC.jl package.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"using SmoQySynthAC\nusing Distributions\nusing CairoMakie\n\n# render figures using SVG backend so they are sharp\nCairoMakie.activate!(type = \"svg\")","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"In this example we will work with the single-particle imaginary time fermion Green's function which is given by","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"G(tau) = int_-infty^infty K(omegataubeta) A(omega)","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"where A(omega) is the spectral function and","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"K(omegataubeta) = frace^-tau omega1 + e^-beta omega","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"is the kernel function where beta = 1T is the inverse temperature and it is assumed that tau in 0 beta).","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"As a first step in demonstrating the functionality of SmoQySynthAC.jl package, let us define a synthetic spectral function A(omega). For convenience we will do this using the Distributions.jl package. We will define a spectral function with a Lorentzian (Cauchy) distribution in centered between two Normal distributions on either side.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"# define spectral distribution\nspectral_dist = MixtureModel(\n    [Normal(-2.0,0.7), Cauchy(0.0, 0.3), Normal(+2.0,0.7)],\n    [0.2, 0.4, 0.4]\n)\n\n# define function to evaluate the spectral funciton\nspectral_function = ω -> pdf(spectral_dist, ω)","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Now let us quickly plot our spectral function so we can see what it looks like.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"ωmin = -6.0\nωmax = 6.0\nω = collect(range(start = ωmin, stop = ωmax, length = 1000))\nA = spectral_function.(ω)\n\nfig = Figure(\n    size = (700, 400),\n    fonts = (; regular= \"CMU Serif\"),\n    figure_padding = 10\n)\n\nax = Axis(\n    fig[1, 1],\n    aspect = 7/4,\n    xlabel = L\"\\omega\", ylabel = L\"A(\\omega)\",\n    xlabelsize = 30, ylabelsize = 30,\n    xticklabelsize = 24, yticklabelsize = 24,\n)\n\nlines!(\n    ω, A,\n    linewidth = 2,\n    alpha = 1.0,\n    color = :black,\n    linestyle = :solid\n)\n\nxlims!(ax, ωmin, ωmax)\nylims!(ax, 0.0, 1.05*maximum(A))\n\nfig","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"The next step is to define the inverse temperature beta, discretization in imaginary Deltatau and corresponding tau grid.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"β = 10.0\nΔτ = 0.05\nτ = collect(range(start = 0.0, stop = β, step = Δτ));\nnothing #hide","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Now we can calculate G(tau) using the spectral_to_imaginary_time_correlation_function method and appropriate kernel funciton kernel_tau_fermi.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Gτ = spectral_to_imaginary_time_correlation_function(\n    τ = τ,\n    β = β,\n    spectral_function = spectral_function,\n    kernel_function = kernel_tau_fermi,\n    tol = 1e-10\n);\nnothing #hide","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Having calculated the exact G(tau) function, let us now add some noise to it using the add_noise method.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Gτ_noisy = add_noise(\n    Cτ_exact = Gτ,\n    τ = τ,\n    σ = 1.0e-2,\n    ξ = 0.5,\n    sum_rule = C0 -> 1 - C0, # enforces fermionic boundary condition that C(τ=β) = 1 - C(τ=0)\n    noise_type = :TruncatedNormal\n);\nnothing #hide","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Let us now plot what both G(tau) and G_rm noisy(tau).","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"fig = Figure(\n    size = (700, 500),\n    fonts = (; regular= \"CMU Serif\"),\n    figure_padding = 10\n)\n\nax = Axis(fig[1, 1],\n    aspect = 7/5,\n    xlabel = L\"\\tau\",\n    ylabel = L\"G_\\sigma(\\tau)\",\n    xlabelsize = 30, ylabelsize = 30,\n    xticklabelsize = 24, yticklabelsize = 24,\n)\n\nxlims!(ax, 0.0, β)\nylims!(ax, 0.0, 1.0)\n\nlines!(\n    τ, Gτ,\n    linewidth = 2, alpha = 2.0, color = :black, linestyle = :solid, label = \"Exact\"\n)\n\nlines!(\n    τ, Gτ_noisy,\n    linewidth = 2, alpha = 1.5, color = :red, linestyle = :solid, label = \"Noisy\"\n)\n\naxislegend(\n    ax, halign = :left, valign = :top, labelsize = 30\n)\n\nfig","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = SmoQySynthAC","category":"page"},{"location":"#SmoQySynthAC","page":"Home","title":"SmoQySynthAC","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"SmoQySynthAC.jl is a package developed to generate synthetic imaginary time correlation data mean to resemble noisy quantum Monte Carlo (QMC) measurements, where the introduced noise is correlated in imaginary time. This synthetic data can then be used to benchmark analytic continuation (AC) methods, or numerical methods for transforming noisy imaginary time correlation data to Matsubara space. The approach used to generate the noisy synthetic data is based method introduced in Shao et. al., the citation for which is given below.","category":"page"},{"location":"","page":"Home","title":"Home","text":"@article{PhysRevX.7.041072,\n  title = {Nearly Deconfined Spinon Excitations in the Square-Lattice Spin-$1/2$ Heisenberg Antiferromagnet},\n  author = {Shao, Hui and Qin, Yan Qi and Capponi, Sylvain and Chesi, Stefano and Meng, Zi Yang and Sandvik, Anders W.},\n  journal = {Phys. Rev. X},\n  volume = {7},\n  issue = {4},\n  pages = {041072},\n  numpages = {26},\n  year = {2017},\n  month = {Dec},\n  publisher = {American Physical Society},\n  doi = {10.1103/PhysRevX.7.041072},\n  url = {https://link.aps.org/doi/10.1103/PhysRevX.7.041072}\n}","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install SmoQySynthAC.jl, simply open the Julia REPL and run the commands","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> ]\npkg> add SmoQySynthAC","category":"page"},{"location":"","page":"Home","title":"Home","text":"or equivalently via Pkg do","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg; Pkg.add(\"SmoQySynthAC\")","category":"page"},{"location":"#Documentation","page":"Home","title":"Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"STABLE: Documentation for the latest version of the code published to the Julia General registry.\nDEV: Documentation for the latest commit to the main branch.","category":"page"},{"location":"#Funding","page":"Home","title":"Funding","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The development of this code was supported by the U.S. Department of Energy, Office of Science, Basic Energy Sciences, under Award Number DE-SC0022311.","category":"page"}]
}
