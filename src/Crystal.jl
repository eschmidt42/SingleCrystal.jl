module Crystal

using LinearAlgebra
using JSON

# https://web.archive.org/web/20080324193801/http://cst-www.nrl.navy.mil/lattice/spcgrp/
# coordinate system -> standard cartesian
# primitive vectors -> describe the primitive cell
# basis vectors -> describe positions of atoms in the primitive cell
# translation of basis vectors by primitive vectors generate a crystal

# https://wiki.fysik.dtu.dk/ase/ase/spacegroup/spacegroup.html

"""
    get_chemical_info()

Returns instances `chemical_symbols`, `atomic_numbers` and `masses` for the setup of `Cell`.
"""
function get_chemical_info()
    chemical_symbols = [
        # 0
        "X",
        # 1
        "X", "He",
        # 2
        "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        # 3
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
        # 4
        "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr",
        # 5
        "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
        "In", "Sn", "Sb", "Te", "I", "Xe",
        # 6
        "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
        "Ho", "Er", "Tm", "Yb", "Lu",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
        "Po", "At", "Rn",
        # 7
        "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
        "Cf", "Es", "Fm", "Md", "No", "Lr",
        "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc",
        "Lv", "Ts", "Og"
    ]

    atomic_numbers = Dict()
    for (Z, symbol) in enumerate(chemical_symbols)
        atomic_numbers[symbol] = Z
    end

    atomic_masses_iupac2016 = [
        1.0,  # X
        1.008,  # H [1.00784, 1.00811]
        4.002602,  # He
        6.94,  # Li [6.938, 6.997]
        9.0121831,  # Be
        10.81,  # B [10.806, 10.821]
        12.011,  # C [12.0096, 12.0116]
        14.007,  # N [14.00643, 14.00728]
        15.999,  # O [15.99903, 15.99977]
        18.998403163,  # F
        20.1797,  # Ne
        22.98976928,  # Na
        24.305,  # Mg [24.304, 24.307]
        26.9815385,  # Al
        28.085,  # Si [28.084, 28.086]
        30.973761998,  # P
        32.06,  # S [32.059, 32.076]
        35.45,  # Cl [35.446, 35.457]
        39.948,  # Ar
        39.0983,  # K
        40.078,  # Ca
        44.955908,  # Sc
        47.867,  # Ti
        50.9415,  # V
        51.9961,  # Cr
        54.938044,  # Mn
        55.845,  # Fe
        58.933194,  # Co
        58.6934,  # Ni
        63.546,  # Cu
        65.38,  # Zn
        69.723,  # Ga
        72.630,  # Ge
        74.921595,  # As
        78.971,  # Se
        79.904,  # Br [79.901, 79.907]
        83.798,  # Kr
        85.4678,  # Rb
        87.62,  # Sr
        88.90584,  # Y
        91.224,  # Zr
        92.90637,  # Nb
        95.95,  # Mo
        97.90721,  # 98Tc
        101.07,  # Ru
        102.90550,  # Rh
        106.42,  # Pd
        107.8682,  # Ag
        112.414,  # Cd
        114.818,  # In
        118.710,  # Sn
        121.760,  # Sb
        127.60,  # Te
        126.90447,  # I
        131.293,  # Xe
        132.90545196,  # Cs
        137.327,  # Ba
        138.90547,  # La
        140.116,  # Ce
        140.90766,  # Pr
        144.242,  # Nd
        144.91276,  # 145Pm
        150.36,  # Sm
        151.964,  # Eu
        157.25,  # Gd
        158.92535,  # Tb
        162.500,  # Dy
        164.93033,  # Ho
        167.259,  # Er
        168.93422,  # Tm
        173.054,  # Yb
        174.9668,  # Lu
        178.49,  # Hf
        180.94788,  # Ta
        183.84,  # W
        186.207,  # Re
        190.23,  # Os
        192.217,  # Ir
        195.084,  # Pt
        196.966569,  # Au
        200.592,  # Hg
        204.38,  # Tl [204.382, 204.385]
        207.2,  # Pb
        208.98040,  # Bi
        208.98243,  # 209Po
        209.98715,  # 210At
        222.01758,  # 222Rn
        223.01974,  # 223Fr
        226.02541,  # 226Ra
        227.02775,  # 227Ac
        232.0377,  # Th
        231.03588,  # Pa
        238.02891,  # U
        237.04817,  # 237Np
        244.06421,  # 244Pu
        243.06138,  # 243Am
        247.07035,  # 247Cm
        247.07031,  # 247Bk
        251.07959,  # 251Cf
        252.0830,  # 252Es
        257.09511,  # 257Fm
        258.09843,  # 258Md
        259.1010,  # 259No
        262.110,  # 262Lr
        267.122,  # 267Rf
        268.126,  # 268Db
        271.134,  # 271Sg
        270.133,  # 270Bh
        269.1338,  # 269Hs
        278.156,  # 278Mt
        281.165,  # 281Ds
        281.166,  # 281Rg
        285.177,  # 285Cn
        286.182,  # 286Nh
        289.190,  # 289Fl
        289.194,  # 289Mc
        293.204,  # 293Lv
        293.208,  # 293Ts
        294.214,  # 294Og
    ]

    masses = Dict(k => atomic_masses_iupac2016[atomic_numbers[k]] for k in chemical_symbols)
    
    return chemical_symbols, atomic_numbers, masses
end

function load_spgs(;file::String="$(@__DIR__)/spacegroup.json")
    return JSON.parsefile(file)
end

function load_refs(;file::String="$(@__DIR__)/ase-atoms.json")
    return JSON.parsefile(file)
end


"""
    parse_json_spacegroup(spgs,nr,setting)

Formats the info of a specific spacegroup, preparing its use in `get_symops`.
"""
function parse_json_spacegroup(spgs::Dict,nr::T0,setting::T1) where {T0<:Integer,T1<:Integer}
    spg = spgs["$(nr): $(setting)"]
    spg["subtrans"] = [Array{Float64}(reshape(v,3)) for v in spg["subtrans"]]
    spg["translations"] = [Array{Float64}(reshape(v,3)) for v in spg["translations"]]
    spg["rotations"] = [Array{Float64}(hcat(v...)) for v in spg["rotations"]]
    return spg
end

"""
    get_symops(spg)

Generates rotation and translation operations for the given spacegroup info (`spg` from `parse_json_spacegroup`).
"""
function get_symops(spg::Dict)
    #  `spg` is the output of `parse_json_spacegroup`
    parities = spg["centrosymmetric"] ? [1,-1] : [1]
    symops = []
    @assert length(spg["rotations"]) == length(spg["translations"])
    for (parity, trans_sub) in Iterators.product(parities, spg["subtrans"])
        for (rot, trans) in Iterators.zip(spg["rotations"], spg["translations"])
            push!(symops, (parity * rot,
                           (trans + trans_sub) .% 1))
        end
    end
    return symops
end

"""
    fold(x)

Folding `x` back into [0,1).
"""
function fold(x::T) where T<:Any
    return x < 0 ? x + 1 : x
end

"""
    get_equivalent_sites(basis, symops)

Computes equivalent sites of in `basis` using rotations and translations in `symops`.
"""
function get_equivalent_sites(basis::Array{Array{Float64,2},1}, symops)
    kinds, sites = [], []

    for (kind, pos) in enumerate(basis)
        for (rot, trans) in symops
            site = (transpose(pos * rot) + trans) .% 1
            site = fold.(site)
            isdifferent = !any([v ≈ site for v in sites])
            if ((length(sites) == 0) | isdifferent)
                push!(sites, site)
                append!(kinds, kind)
            end
        end
    end
    return sites, kinds
end

"""
    make_unit_vec(x)

Computes the normalized `x` vector.
"""
function make_unit_vec(x::Array{T,1}) where T<:Any
    return x / norm(x)
end

"""
    get_coords(a_direction, ab_normal)

Computes the coordinate system given `a_direction` (usually [1 0 0]) and the vector normal to `a` and `b` (`ab_normal`, usually [0 0 1]).
"""
function get_coords(a_direction::Array{Float64,1}, ab_normal::Array{Float64,1})
    @assert dot(a_direction, ab_normal) ≈ 0.
    _x = make_unit_vec(a_direction)
    z = make_unit_vec(ab_normal)

    x = _x - dot(_x, ab_normal) * z
    xyz = hcat(x, cross(z,x), z)
    return xyz
end

"""
    deg2rad(x)

Convert degrees to radians.
"""
function deg2rad(x::T) where T <: Real
    return x * π / 180.
end

"""
    get_cos(x)

Compute the cosine of x degrees.
"""
function get_cos(x::T) where T <: Real
    return x ≈ 90 ? 0 : cos(deg2rad(x))
end

"""
    get_cell_vectors(cellpar)

Computes 3x3 cell vectors (column-wise). cellpar = [a, b, c, α, β, γ]
"""
function get_cell_vectors(cellpar::Array{Float64,1})
    a, b, c, α, β, γ = cellpar
    cos_α = get_cos(α)
    cos_β = get_cos(β)
    cos_γ = get_cos(abs(γ))
    sin_γ = abs(γ) ≈ 90 ? sign(γ) : sin(deg2rad(γ))
    cos_α, cos_β, cos_γ, sin_γ
    
    cy = (cos_α - cos_β * cos_γ) / sin_γ
    abc = hcat([a; 0; 0], b*[cos_γ; sin_γ; 0], c*[cos_β; cy; √(1-cos_β*cos_β-cy*cy)])
    return abc
end


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
    Atom(; name, mass)

Generates an atom object for the creation of unit cells.
"""
struct Atom{T}
    name::String
    mass::T
end

Atom(name,mass) = Atom{Float64}(name,mass)
make_atom(name,mass) = Atom{typeof(mass)}(name,mass)

function Atom(;name::String="dummy",mass::T=42) where T<:Real
    return Atom{typeof(mass)}(name, mass)
end

struct Spacegroup
    nr::Int32 #spacegroup number
    setting::Int8 # spacegroup setting (1 or 2, inherited from ase.spacegroup)
    kinds::Array{Int16,1} 
    sites::Array{Any,1}
end

Spacegroup() = Spacegroup(-1,-1,[-1],[0])

struct Cell{T}
    atoms::Array
    positions
    box::CoordinateSystem{T}
    edge_lengths::Array{T}
    # spacegroup related infos
    spacegroup::Spacegroup
end

Cell(atoms, positions, box, edge_lengths) = Cell{Float64}(atoms, positions, box, edge_lengths, Spacegroup())
Cell(atoms, positions, box, edge_lengths, spg) = Cell{Float64}(atoms, positions, box, edge_lengths, spg)

"""
    make_unitcell(basis, symbols, nr, setting, cellpar; spg_file=nothing, a_direction=[1.; 0.; 0.], ab_normal=[0.;0.;1.], atom_cb=make_atom)

Creates a `Cell` struct including atom information, atom coordinates, the cell vectors/box and the edge lengths of the box. 
"""
function make_unitcell(basis, symbols, nr, setting, cellpar; 
        a_direction=[1.; 0.; 0.], ab_normal=[0.; 0.; 1.], 
        atom_cb::Function=make_atom)
    # get constants
    chemical_symbols, atomic_numbers, masses = get_chemical_info() 

    # get stored symmetry operation components
    spgs = load_spgs()
    spg = parse_json_spacegroup(spgs, nr, setting)
    
    # prepare symmetry operations for application to `basis`
    symops = get_symops(spg)

    # generating equivalent sites
    sites, kinds = get_equivalent_sites(basis, symops)

    # computing the simulation cell/box
    xyz = get_coords(a_direction, ab_normal)
    abc = get_cell_vectors(cellpar)
    cell = abc * xyz

    # storing crystal properties in a `Crystal.Cell` instance
    el2atom_map = Dict(el => atom_cb(el, masses[el]) for el in symbols)

    cc = CartesianCoords(Float64)
    box = PrimitiveVectors(cc, A₁=abc[:,1], A₂=abc[:,2], A₃=abc[:,3])

    _spg = Spacegroup(nr, setting, kinds, sites)

    crystal = Cell(
        [el2atom_map[symbols[v]] for v in kinds],
        [cell * v for v in  sites],
        box,
        [norm(abc[:,1]), norm(abc[:,2]), norm(abc[:,3])],
        _spg
    )
    return crystal
end

"""
    make_bcc_unitcell(symbol, a)

Convenience function wrapping `make_unitcell` to create a bcc unit cell.
"""
make_bcc_unitcell(symbol::String,a::T,atom_cb::Function=make_atom) where T<:Real = make_unitcell([[0. 0. 0.],], [symbol], 229, 1, [a,a,a,90.,90,90], atom_cb=atom_cb)

"""
    make_fcc_unitcell(element, a)

Convenience function wrapping `make_unitcell` to create a fcc unit cell.
"""
make_fcc_unitcell(symbol::String,a::T,atom_cb::Function=make_atom) where T<:Real = make_unitcell([[0. 0. 0.],], [symbol], 225, 1, [a,a,a,90.,90,90], atom_cb=atom_cb)


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
    positions_vac = [c.positions[i] for i in 1:n_atoms if !(i in ixs)]
    type = typeof(positions_vac[1][1])
    return Cell(atoms_vac, positions_vac, c.box, c.edge_lengths, c.spacegroup)
end

"""
    make_supercell(c;nx=1,ny=1,nz=1)

Creates a supercell by replicating the unitcell `c` and translating it along the vectors in `cell.box` `nx` * `ny` * `nz` times.
"""
function make_supercell(c::Cell; nx::I=1, ny::I=1, nz::I=1) where {I <: Integer}
    @assert (nx > 0) & (ny > 0) & (nz > 0) 
    atoms = []
    positions = []
    box = PrimitiveVectors(c.box, A₁=[nx; 0; 0], A₂=[0; ny; 0], A₃=[0; 0; nz])
    edge_lengths = [norm(box.x), norm(box.y), norm(box.z)]

    for i in 0:nx-1, j in 0:ny-1, k in 0:nz-1
        shift = c.box.M * [i, j, k] # [c.box.x * i + c.box.y * j + c.box.z * k]
        new_pos = [pos + shift for pos in c.positions]
        push!(positions, new_pos) # TODO:  ((3, 1), Size(3,)) of input arrays do not match
        push!(atoms, c.atoms)
    end
    supercell = Cell(vcat(atoms...), vcat(positions...), box, edge_lengths)
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

function parse_json_crystal(json_ase_crystal::Dict{String,Any})
    ase_crystal = Dict()
    ase_crystal["spacegroup_kinds"] = Array{Int32}(json_ase_crystal["spacegroup_kinds"])
    ase_crystal["positions"] = [Array{Float64}(v) for v in json_ase_crystal["positions"]]
    ase_crystal["cell"] = hcat(json_ase_crystal["cell"]...)
    ase_crystal["symbols"] = Array{String}(json_ase_crystal["symbols"])
    return ase_crystal
end

function cell2dict(crystal::Crystal.Cell)

    chemical_symbols, atomic_numbers, masses = get_chemical_info()

    d = Dict() 
    d["numbers"] = [atomic_numbers[atom.name] for atom in crystal.atoms] #[atomic_numbers[symbols[v]] for v in kinds]
    d["positions"] = hcat(crystal.positions...)
    d["spacegroup_kinds"] = crystal.spacegroup.kinds
    d["cell"] = crystal.box.M
    d["pbc"] = [true,true,true]
    d["info"] = Dict("spacegroup"=>Dict("nr"=>crystal.spacegroup.nr, "setting"=>crystal.spacegroup.setting))
    return d
end

end