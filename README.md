# Representation of Magnetic Fields using Solid Harmonic Expansions - Example

This repository contains example code to represent a magnetic field with a solid spherical harmonic expansion. Therefore, a spherical t-design is measured on a sphere to calculate the coefficients of the expansion by an efficient quadrature. Then, the magnetic field can be described by 
```math 
        \boldsymbol B^{\boldsymbol\rho}(\boldsymbol a) = 
        \sum_{l=0}^L\sum_{m=-l}^l \boldsymbol\gamma_{l,m}(\boldsymbol\rho) Z_l^m(\boldsymbol a)
        \qquad \forall \boldsymbol a\in\mathcal{B}_R(\boldsymbol 0), 
```
where $\boldsymbol B^{\boldsymbol\rho}:\Omega \subseteq\mathbb{R}^3 \rightarrow \mathbb{R}^3$ is the magnetic field with indicated origin of the coordinate system $\boldsymbol\rho\in\Omega$, $Z_l^m : \mathbb{R} \rightarrow \mathbb{R}$ are the normalized real solid spherical harmonics, and $\mathcal B_R(\boldsymbol 0)$ is a ball with radius $R$ defined the measurement. Finally, $\boldsymbol\gamma_{l,m}(\boldsymbol\rho)\in\mathbb{R}^3$ are the coefficients, which are calculated in this code example.


The theory and methods are described in the associated publication

M. Boberg, T. Knopp, and M. Möddel. *(under review)*. Unique Compact Representation of Magnetic Fields using Truncated Solid Harmonic Expansions. 



## Installation

In order to use this code one first has to download [Julia](https://julialang.org/) (version 1.8 or later), clone this repository and navigate to the folder in the command line. The example script automatically activates the environment and installs all necessary packages.

## Execution
After installation the example code can be executed by running `julia` in the REPL and entering
```julia
include("example.jl")
```
The example is structured into three parts:
1. Loading the measured data of a magnetic gradient field at $36$ spherical 8-design positions and calculating the initial coefficients by an efficient quadrature.
2. Post processing the coefficients: 
    * Correcting the offset in the Hall-effect sensor by shifting the coefficients into a common coordinate system.
    * Calculating the field-free-point of the gradient field and shifting the coefficients into this unique point, i.e., changing the origin of the coordinate system $\boldsymbol\rho$.
    * Additionally, an error propagation of the measurement errors is provided.
3. Visualizing of the final coefficients and magnetic field.

The main functions for calculating and shifting the coefficients are provided by the package [SphericalHarmonicExpansions.jl](https://github.com/hofmannmartin/SphericalHarmonicExpansions.jl). All functions that are magnetic field specific (like the calculation of the field-free-point or visualization) are provided in the utils folder.


## Open MPI Data

The measurement data associated to this project stems from a $2~\textup{T m}^{-1}$ magnetic gradient field used for spatial encoding in magnetic particle imaging (MPI). The data is stored in the [HDF5 file format](https://www.hdfgroup.org/solutions/hdf5/) and can be found in the `data` folder.
