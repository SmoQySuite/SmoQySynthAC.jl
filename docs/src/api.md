# API

## Spectral to Correlation Function

- [`spectral_to_imaginary_time_correlation_function`](@ref)
- [`spectral_to_matsubara_correlation_function`](@ref)

```@docs
spectral_to_imaginary_time_correlation_function
spectral_to_matsubara_correlation_function
```

## Add Noise to an Imaginary Time Correlation Function

- [`add_noise`](@ref)
- [`add_noise!`](@ref)

```@docs
add_noise
add_noise!
```

## Kernel Functions

- [`kernel_mat`](@ref)
- [`kernel_tau_fermi`](@ref)
- [`kernel_mat_fermi`](@ref)
- [`kernel_tau_bose`](@ref)
- [`kernel_mat_bose`](@ref)
- [`kernel_tau_bose_alt`](@ref)
- [`kernel_mat_bose_alt`](@ref)
- [`kernel_tau_sym_bose`](@ref)
- [`kernel_mat_sym_bose`](@ref)
- [`kernel_tau_sym_bose_alt`](@ref)
- [`kernel_mat_sym_bose_alt`](@ref)

```@docs
kernel_mat
```

### Fermionic Kernel Functions

```@docs
kernel_tau_fermi
kernel_mat_fermi
```

### Bosonic Kernel Functions

```@docs
kernel_tau_bose
kernel_mat_bose
```

### Modified Bosonic Kernel Functions

```@docs
kernel_tau_bose_alt
kernel_mat_bose_alt
```

### Symmetric Bosonic Kernel Function

```@docs
kernel_tau_sym_bose
kernel_mat_sym_bose
```

### Modified Symmetric Bosonic Kernel Function

```@docs
kernel_tau_sym_bose_alt
kernel_mat_sym_bose_alt
```

## Quantum Statistics Functions

- [`fermi`](@ref)
- [`bose`](@ref)

```@docs
fermi
bose
```