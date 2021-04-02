# Crystal.jl

This package contains functionality to create body centered cubic and face centered cubic crystal, both of which are very common structures in metals ([💡 wiki 💡](https://en.wikipedia.org/wiki/Cubic_crystal_system)).

When writing this my main objective was to create input for **[Molly.jl](https://github.com/JuliaMolSim/Molly.jl)**, a molecular dynamics simulations package. If you don't know about it yet, I highly recommend to check it out! 😁 

## Installation

To install, open your Julia REPL and enter
```Julia
import Pkg
Pkg.add("https://github.com/eschmidt42/Crystal.jl")
```

or alternatively via the package manager (enter in the REPL via `]`) and type `add https://github.com/eschmidt42/Crystal.jl`.

## Usage

`docs/crystal.ipynb` describes how to use this package, going from the creation of a unit cell to the 3d animation of a supercell with a vacancy (missing atom).

To run this notebook you'll need the acompanying `docs/Project.toml` and `docs/Manifest.toml` files. In case you download this package just `cd` into `docs` using Jupyter and execute the notebook one cell at a time. 

You may notice that the notebook uses additional packages, which are not required for functions of `Crystal.jl` itself, like `DataFrames`. Those packages help to visualize the results.   
