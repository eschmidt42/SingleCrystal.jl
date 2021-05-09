# SingleCrystal.jl

With this package you can create single crystal structures. The implementation is based on the Python [**ase**](https://gitlab.com/ase/ase) package.

## Installation

To install, open your Julia REPL and enter
```Julia
import Pkg
Pkg.add("https://github.com/eschmidt42/SingleCrystal.jl#master")
```

or alternatively via the package manager (enter in the REPL via `]`) and type `add https://github.com/eschmidt42/SingleCrystal.jl#master`.

## Usage

Let's say you want to create a body centered cubic (bcc) unit cell ([💡 wiki 💡](https://en.wikipedia.org/wiki/Cubic_crystal_system)). The general approach to create it or any other crystal's unit cell would be as follows:

```Julia
symbols = ["Fe"] # chemical elements
basis = [[0. 0. 0.],] # scaled coordinates
nr = 229 # space group
setting = 1 # space group settig (greetings from ase)
cellpar = [2.87, 2.87, 2.87, 90, 90, 90]; # specification of the 3 cell vector lengths (in Å = 10⁻¹⁰m) a, b, c and three angles (in degrees) α, β, γ

crystal = SingleCrystal.make_unitcell(basis, symbols, nr, setting, cellpar)
```

For the case of bcc unit cells, you could alternatively also use the less verbose path via the `make_bcc_unitcell` convenience function:

```Julia
crystal = SingleCrystal.make_bcc_unitcell("Fe", 3.4)
```

In case you want to replicate the unit cell along the cell vectors to create a supercell, you can use `make_supercell`:

```Julia
supercell = SingleCrystal.make_supercell(crystal, nx=3, ny=3, nz=3);
```

For more examples, and a peek behind the curtains of the ase algorithm implemented in this package, I encourage you to check out `docs/singl_crystals_in_julia.ipynb`. There you can also find the above examples in context and find how to create a vacancy \*spoiler\*.

Happy crystal synthesizing! 😃

## Motivation

The main objective for this package is to prepare input required for the Molecular Dynamics package **[Molly.jl](https://github.com/JuliaMolSim/Molly.jl)**, to simulate body centered cubic single crystals. 

A minimal working example for the usage of SingleCrystal.jl with Molly.jl (based on a fork of Molly.jl adding the Finnis-Sinclair potential type - a pull request of the fork is currently [under review](https://github.com/JuliaMolSim/Molly.jl/pull/32)):

```Julia
import Pkg
Pkg.activate(".") # if you are in the root of the forked Molly.jl package

using Molly
using SingleCrystal

fs_inter, elements, masses, bcc_lattice_constants, reference_energies = Molly.get_finnissinclair1984(true)
make_atom(name,mass) = Atom(name=name,mass=mass)

# setting up the crystal
nx = 3
ny = 3
nz = 3
element = "Fe"

a = bcc_lattice_constants[element]
crystal = SingleCrystal.make_bcc_unitcell(element, a, make_atom)
supercell = SingleCrystal.make_supercell(crystal, nx=nx, ny=ny, nz=nz)

# setting up the simulation
T = 100. # Kelvin
T = T*fs_inter.kb
n_steps = 2500
dt = .002 # ps; ns = 1e-9s, ps = 1e-12s, fs = 1e-15s
n_atoms = length(supercell.atoms)
general_inters = (fs_inter,)
velocities = [velocity(supercell.atoms[i].mass, T, dims=3) for i in 1:n_atoms]
nb_matrix = trues(n_atoms,n_atoms)
dist_cutoff = 2 * a
nf = DistanceNeighbourFinder(nb_matrix, 1, dist_cutoff)
thermostat = NoThermostat()

loggers = Dict(
    "temperature" => TemperatureLogger(1),
    "pot" => EnergyLogger(1),
)

s = Simulation(
    simulator=VelocityVerlet(), 
    atoms=supercell.atoms, 
    general_inters=general_inters,
    coords=[SVector{3}(v) for v in supercell.positions], 
    velocities=velocities,
    temperature=T, 
    box_size=supercell.edge_lengths[1],
    timestep=dt,
    n_steps=n_steps,
    neighbour_finder=nf,
    loggers=loggers,
)

# running the simulation
simulate!(s) 
```

## To dos

1. Add more crystal structures to test beyond those in the ase gallery
2. Test how the performance / resource requirements scale with crystal size / number of atoms

## Contributing

Contributions are very welcome. 