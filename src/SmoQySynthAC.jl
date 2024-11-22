module SmoQySynthAC

using FFTW
using QuadGK
using Optim
using Distributions
using SpecialFunctions

# generate noise based on a truncated normal distribution define on the domain (0, ∞)
include("truncated_normal_noise.jl")

# define various types of kernel functions used for analytic continuation
include("kernel_functions.jl")
export kernel_tau_fermi, kernel_tau_bose, kernel_tau_sym_bose
export kernel_mat_fermi, kernel_mat_bose
export fermi, bose

# functions that take as inputs the spectral function and return various correlation kernel_functions
include("spectral_to_correlation_function.jl")
export spectral_to_imaginary_time_correlation_function
export spectral_to_matsubara_correlation_function

# functions for adding noise to imaginary time correlation function
include("add_noise.jl")
export add_noise!, add_noise

end
