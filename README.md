# CARP
Chandrasekhar Anisotropic Rotating Plummer.

This repository allows for the computation of the diffusion coefficients and the diffusion rates of a rotating, anisotropic Plummer cluster using the orbit-averaged Chandrasekhar theory.

The parametrization of velocity anisotropy follows [Dejonghe (1987)](https://ui.adsabs.harvard.edu/abs/1987MNRAS.224...13D/abstract), and is described by a parameter $q$.

The parametrization rotation follows the Lynden-Bell demon from [Lynden-Bell (1960)](https://ui.adsabs.harvard.edu/abs/1960MNRAS.120..204L/abstract), and described by

$$F_{\mathrm{rot}}(J_{\mathrm{r}},L,L_{z}) = F_{\mathrm{tot}}(J_{\mathrm{r}},L) \big(1 + \alpha \, \mathrm{sgn}[L_{\mathrm{z}}/L] \big),$$

where $F_{\mathrm{rot}}(J_{\mathrm{r}},L,L_{z})$ is the distribution function in action space for the rotating cluster, $F_{\mathrm{tot}}(J_{\mathrm{r}},L)$ is the distribution function for the associated non-rotating cluster and $\alpha$ is the rotation parameter (between -1 and 1).

This repository makes use of the Julia Language  (see [Bezanson et al. (2017)](https://doi.org/10.1137/141000671)), whose repository is located at [JuliaLang](https://github.com/JuliaLang/julia/tree/master).

## Installation

Install Julia by following the instructions at [julialang.org/downloads/](https://julialang.org/downloads/).

To invoke Julia in the Terminal, you need to make sure that the julia command-line program is in your `PATH`. See [here](https://julialang.org/downloads/platform/#optional_add_julia_to_path) for detailed instructions.

## Packages

Open the terminal in the folder `packages` and type

```
$ julia Install-pkg.jl
```

to install the following packages:

- `ArgParse`
- `HypergeometricFunctions`
- `SpecialFunctions`
- `StaticArrays`
- `Interpolations`
- `HDF5`

### !! WARNING !!

**DO NOT INTERRUPT THE DOWNLOADING OF THE PACKAGES !!!!**


## Sampling the diffusion rate in action space

To compute a mapping of the diffusion rate $\partial F/\partial t$, one can use the Julia scripts located at
`code/compute/`. These script will compute a 2D map in action space and store it in a `.hf5` file in the folder `code/data`. 

As an example, `Compute_dFdt_JrL.jl` will compute a 2D map in $(J_{\mathrm{r}},L)$ space.

This script is parallelized. The full list of arguments can be found in the file `code/sources/julia/Args.jl`.

Here is an example of a console command to launch the computation in $(J_{\mathrm{r}},L)$ space for 12 threads , for a Plummer cluster with anisotropy parameter $q=0.0$ and rotation parameter $\alpha=0.1$.

```
$ cd ./code/compute
$ julia -t 12 Compute_dFdt_JrL.jl --q 0.0 --a 0.1
```


The resulting file will be created in the folder `code/data` under the name 
`Dump_dFdt_JrL_Map_q_0.0_a_0.1.hf5`.
