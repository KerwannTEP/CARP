# CARP
Chandrasekhar Anisotropic Rotating Plummer.

## Aim of the code

This repository allows for the computation of the diffusion coefficients and the diffusion rates of a rotating, anisotropic Plummer cluster using the orbit-averaged Chandrasekhar theory.

The parametrization of velocity anisotropy follows [Dejonghe (1987)](https://ui.adsabs.harvard.edu/abs/1987MNRAS.224...13D/abstract), and is described by a parameter $q$.

The parametrization rotation follows the Lynden-Bell demon from [Lynden-Bell (1960)](https://ui.adsabs.harvard.edu/abs/1960MNRAS.120..204L/abstract), and described by

$$F_{\mathrm{rot}}(J_{\mathrm{r}},L,L_{z}) = F_{\mathrm{tot}}(J_{\mathrm{r}},L) \big(1 + \alpha \ \mathrm{sgn}[L_{\mathrm{z}}/L] \big),$$

where $F_{\mathrm{rot}}(J_{\mathrm{r}},L,L_{z})$ is the distribution function in action space for the rotating cluster, $F_{\mathrm{tot}}(J_{\mathrm{r}},L)$ is the distribution function for the associated non-rotating cluster and $\alpha$ is the rotation parameter (between -1 and 1).

This repository makes use of the Julia Language  (see [Bezanson et al. (2017)](https://doi.org/10.1137/141000671)), whose repository is located at [JuliaLang](https://github.com/JuliaLang/julia/tree/master).

A GPU-accelerated version of this code, making use of the `CUDA.jl` library, is available at the repository [CARP_GPU](https://github.com/KerwannTEP/CARP_GPU).

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
- `Plots`
- `LaTeXStrings`

### !! WARNING !!

**DO NOT INTERRUPT THE DOWNLOADING OF THE PACKAGES !!!!**


## Sampling the diffusion rate in action space

To compute a mapping of the diffusion rate $\partial F/\partial t$, one can use the Julia scripts located at
`code/compute/`. These script will compute a 2D map in action space and store it in a `.hf5` file in the folder `code/data`. 

As an example, `Compute_dFdt_JrL.jl` will compute a 2D map in $(J_{\mathrm{r}},L)$ space.

This script is parallelized. The full list of arguments can be found in the file `code/sources/julia/Args.jl`.

Here is an example of a console command to launch the computation in $(J_{\mathrm{r}},L)$ space for 12 threads, for a Plummer cluster with anisotropy parameter $q=0.0$ and rotation parameter $\alpha=0.1$.

```
$ cd ./code/compute
$ julia -t 12 Compute_dFdt_JrL.jl --q 0.0 --a 0.1
```


The resulting file will be created in the folder `code/data` under the name 
`Dump_dFdt_JrL_Map_q_0.0_a_0.1.hf5` in about 30 seconds for these parameters.


One may also be interested in samlping the $(J_{\mathrm{r}},\cos I)$ space. To that end, one should use the following command

```
$ cd ./code/compute
$ julia -t 12 Compute_dFdt_JrCosI.jl --q 0.0 --a 0.1 --nbw 20 --nbphi 200 --nbvartheta 50
```

The resulting file will be created in the folder `code/data` under the name 
`Dump_dFdt_JrCosI_Map_q_0.0_a_0.1.hf5` in about 120 seconds for these parameters.

## Plotting the diffusion rate in action space

After having sampled a 2D-slice of the action space, one may plot the data using the Julia scripts located at `code/compute/`.

As an example, `PlotContour_JrL.jl` will compute a 2D map in $(J_{\mathrm{r}},L)$ space for the data of the cluster $(q,\alpha)=(0.0,0.1)$ when typing the follow command in the console

```
$ cd ./code/compute
$ julia PlotContour_JrL.jl --q 0.0 --a 0.1
```

The resulting file will be created in the folder `code/graphs` under the name 
`Map_dFdt_JrL_q_0.0_a_0.1.png`.


![`dF/dt (Jr,L)` for `q=0` and `$alpha=0.1`](code/graphs/examples/Map_dFdt_JrL_q_0.0_a_0.1.png)

Similarly, `PlotContour_JrCosI.jl` will compute a 2D map in $(J_{\mathrm{r}},\cos I)$ space for the data of the cluster $(q,\alpha)=(0.0,0.1)$ when typing the follow command in the console

```
$ cd ./code/compute
$ julia PlotContour_JrCosI.jl --q 0.0 --a 0.1
```

The resulting file will be created in the folder `code/graphs` under the name 
`Map_dFdt_JrCosI_q_0.0_a_0.1.png`.


![`dF/dt (Jr,cos I)` for `q=0` and `$alpha=0.1`](code/graphs/examples/Map_dFdt_JrCosI_q_0.0_a_0.1.png)
