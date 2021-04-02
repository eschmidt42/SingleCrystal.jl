module Crystal

using LinearAlgebra

# https://web.archive.org/web/20080324193801/http://cst-www.nrl.navy.mil/lattice/spcgrp/
# coordinate system -> standard cartesian
# primitive vectors -> describe the primitive cell
# basis vectors -> describe positions of atoms in the primitive cell
# translation of basis vectors by primitive vectors generate a crystal

# https://wiki.fysik.dtu.dk/ase/ase/spacegroup/spacegroup.html


struct CoordinateSystem{T}
    x::Array{T}
    y::Array{T}
    z::Array{T}
    M::Matrix{T} # contains x,y,z as column vectors
end

"""
    CartesianCoords(T)

Initialization of the cartesian coordiantes using `CoordianteSystem`.
"""
CartesianCoords(T) = CoordinateSystem{T}([1; 0; 0], [0; 1; 0], [0; 0; 1], I(3))

"""
    PrimitiveVectors(cc)

Generates vectors for a primitive cell on the basis of a coordinate system, like the `CartesianCoords`.
"""
function PrimitiveVectors(cc::CoordinateSystem{T}; 
        A₁=[1; 0; 0], A₂=[0; 1; 0], A₃=[0; 0; 1]) where {T<:Real}
    # x = SVector{3}(sum([cc.x, cc.y, cc.z] .* A₁))
    x = cc.M * A₁
    y = cc.M * A₂
    z = cc.M * A₃
    M = hcat(x,y,z)
    return CoordinateSystem{T}(x, y, z, M)
end

"""
    get_basis_vectors(type)

Returns basis vectors for atom locations in a unit cell.
"""
function get_basis_vectors(type::String)
    # type of crystal
    db = Dict(
        "fcc" => [[0 0 0], [0 1//2 1//2], [1//2 0 1//2], [1//2 1//2 0]],
        "bcc" => [[0 0 0], [1//2 1//2 1//2]]
    )
    return db[type]
end

"""
    Atom(; name, mass)

Generates an atom object for the creation of unit cells.
"""
struct Atom{T}
    name::String
    mass::T
end

function Atom(;name::String="dummy",mass::T=42) where T<:Real
    return Atom{typeof(mass)}(name, mass)
end

struct Cell{T}
    atoms::Array
    coords
    box::CoordinateSystem{T}
    edge_lengths::Array{T}
end

"""
    make_unitcell(elements, a, el2atom_map, basis_info)

Creates a `Cell` struct including atom information, atom coordinates, the cell vectors/box and the edge lengths of the box. 
"""
function make_unitcell(elements::Array{String},a::T,el2atom_map::Dict{String},basis_info::String) where T <:Real
    cc = CartesianCoords(T)
    box = PrimitiveVectors(cc, A₁=[a; 0; 0], A₂=[0; a; 0], A₃=[0; 0; a])
    basis = get_basis_vectors(basis_info)
    coords = [box.M * transpose(b) for b in basis]

    @assert length(elements) == length(coords) "This function generates a $(basis_info) unit cell with $(length(coords)) atoms, so 4 strings are required but got $(elements)"
    atoms = [el2atom_map[el] for el in elements]
    edge_lengths = [norm(box.x), norm(box.y), norm(box.z)]
    
    return Cell{T}(atoms, coords, box, edge_lengths)
end

"""
    make_bcc_unitcell(elements, a, el2atom_map)

Convenience function wrapping `make_unitcell` to create a bcc unit cell.
"""
make_bcc_unitcell(elements::Array{String},a::T,el2atom_map::Dict{String}) where T<:Real = make_unitcell(elements,a,el2atom_map,"bcc")

"""
    make_fcc_unitcell(elements, a, el2atom_map)

Convenience function wrapping `make_unitcell` to create a fcc unit cell.
"""
make_fcc_unitcell(elements::Array{String},a::T,el2atom_map::Dict{String}) where T<:Real = make_unitcell(elements,a,el2atom_map,"fcc")

"""
    add_vacancies(c, ixs; random=false, n_vac=1)

Adds a vacancy / removes an atom to/from a `Cell` struct.
"""
function add_vacancies(c::Cell; ixs::Array{Int}=[1], random::Bool=false, n_vac::Int=1)
    # if `random=true` then `n_vac` are generated and `ixs` is ignored
    n_atoms = length(c.atoms)
    if random == true
        ixs = rand(1:n_atoms, n_vac)
    end
    atoms_vac = [c.atoms[i] for i in 1:n_atoms if !(i in ixs)]
    coords_vac = [c.coords[i] for i in 1:n_atoms if !(i in ixs)]
    type = typeof(coords_vac[1][1])
    return Cell{type}(atoms_vac, coords_vac, c.box, c.edge_lengths)
end

"""
    make_supercell(c;nx=1,ny=1,nz=1)

Creates a supercell by replicating the unitcell `c` and translating it along the vectors in `cell.box` `nx` * `ny` * `nz` times.
"""
function make_supercell(c::Cell; nx::I=1, ny::I=1, nz::I=1) where {I <: Integer}
    @assert (nx > 0) & (ny > 0) & (nz > 0) 
    atoms = []
    coords = []
    box = PrimitiveVectors(c.box, A₁=[nx; 0; 0], A₂=[0; ny; 0], A₃=[0; 0; nz])
    edge_lengths = [norm(box.x), norm(box.y), norm(box.z)]

    for i in 0:nx-1, j in 0:ny-1, k in 0:nz-1
        shift = c.box.M * [i, j, k] # [c.box.x * i + c.box.y * j + c.box.z * k]
        new_pos = [pos + shift for pos in c.coords]
        push!(coords, new_pos) # TODO:  ((3, 1), Size(3,)) of input arrays do not match
        push!(atoms, c.atoms)
    end
    supercell = Cell(vcat(atoms...), vcat(coords...), box, edge_lengths)
    return supercell
end

"""
    vector1D(c1,c2,box_size)

Computes the periodic 1D distance between `c1` and `c2` with `box_size` as the periodic boundary.
"""
function vector1D(c1::Real, c2::Real, box_size::Real)
    if c1 < c2
        return (c2 - c1) < (c1 - c2 + box_size) ? (c2 - c1) : (c2 - c1 - box_size)
    else
        return (c1 - c2) < (c2 - c1 + box_size) ? (c2 - c1) : (c2 - c1 + box_size)
    end
end

"""
    vector(c1,c2,box_size)

Convenience function wrapping `vector1D` for the case of a cube.
"""
vector(c1, c2, box_size::Real) = vector1D.(c1, c2, box_size) # all box vectors orthogonal and of equal length

"""
    vector(c1,c2,edge_lengths)

Convenience function wrapping `vector1D` for the case of a non-cubic box with still orthogonal box vectors.
"""
vector(c1, c2, edge_lengths::Array) = vector1D.(c1, c2, edge_lengths) # all box vectors orthogonal

end