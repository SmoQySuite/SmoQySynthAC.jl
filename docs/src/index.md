```@meta
CurrentModule = SmoQySynthAC
```

# SmoQySynthAC

[SmoQySynthAC.jl](https://github.com/SmoQySuite/SmoQySynthAC.jl) is a package developed to generate synthetic imaginary time correlation
data mean to resemble noisy quantum Monte Carlo (QMC) measurements, where the introduced noise is correlated in imaginary time.
This synthetic data can then be used to benchmark analytic continuation (AC) methods, or numerical methods for transforming noisy
imaginary time correlation data to Matsubara space. The approach used to generate the noisy synthetic data is based method introduced
in Shao et. al., the citation for which is given below.

```bibtex
@article{PhysRevX.7.041072,
  title = {Nearly Deconfined Spinon Excitations in the Square-Lattice Spin-$1/2$ Heisenberg Antiferromagnet},
  author = {Shao, Hui and Qin, Yan Qi and Capponi, Sylvain and Chesi, Stefano and Meng, Zi Yang and Sandvik, Anders W.},
  journal = {Phys. Rev. X},
  volume = {7},
  issue = {4},
  pages = {041072},
  numpages = {26},
  year = {2017},
  month = {Dec},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevX.7.041072},
  url = {https://link.aps.org/doi/10.1103/PhysRevX.7.041072}
}
```

## Installation

To install [`SmoQySynthAC.jl`](https://github.com/SmoQySuite/SmoQySynthAC.jl.git),
simply open the Julia REPL and run the commands
```julia
julia> ]
pkg> add SmoQySynthAC
```
or equivalently via `Pkg` do
```julia
julia> using Pkg; Pkg.add("SmoQySynthAC")
```

## Documentation

- [STABLE](https://SmoQySuite.github.io/SmoQySynthAC.jl/stable/): Documentation for the latest version of the code published to the Julia [`General`](https://github.com/JuliaRegistries/General.git) registry.
- [DEV](https://SmoQySuite.github.io/SmoQySynthAC.jl/dev/): Documentation for the latest commit to the `main` branch.

## Funding

The development of this code was supported by the U.S. Department of Energy, Office of Science, Basic Energy Sciences,
under Award Number DE-SC0022311.