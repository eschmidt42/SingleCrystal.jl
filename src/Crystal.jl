module Crystal

using DataFrames
using Plots
using Test
using LaTeXStrings
using LinearAlgebra
using SparseArrays
using Random

masses = Dict("V" => 50.9415, "Nb" => 92.9064, "Ta" => 180.9479,
              "Cr" => 51.996, "Mo" => 95.94, "W" => 183.85,
              "Fe" => 55.847)

struct Atom{T}
    name::String
    mass::T
end

function Atom(;name::String="dummy",mass::T=42) where T<:Real
    return Atom{typeof(mass)}(name, mass)
end

function make_fcc_unitcell(elements::Array{String};a::T=1,el2atom_map::Dict) where T <:Real
    coords = [[0 0 0],[1//2 1//2 0],[1//2 0 1//2],[0 1//2 1//2]]
    @assert length(elements) == length(coords) "This function generates a fcc unit cell with $(length(coords)) atoms, so 4 strings are required but got $(elements)"
    atoms = [el2atom_map[el] for el in elements]
    box_size = Diagonal([a, a, a])
    box_vectors = [1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
    box = box_vectors * box_size
    coords = [v*box for v in coords]
    return atoms, coords, box, box_size, box_vectors
end

function make_bcc_unitcell(elements::Array{String};a::T=1,el2atom_map::Dict) where T <:Real
    coords =[[0 0 0], [.5 .5 .5]]
    @assert length(elements) == length(coords) "This function generates a bcc unit cell with $(length(coords)) atoms, so 4 strings are required but got $(elements)"
    atoms = [el2atom_map[el] for el in elements]
    box_size = Diagonal([a, a, a])
    box_vectors = [1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
    box = box_vectors * box_size
    coords = [v*box for v in coords]
    return atoms, coords, box, box_size, box_vectors
end

function add_vacancies(
        atoms::Array, coords::Array;
        ixs::Array{I}=[2],
        random::Bool=false, 
        n_vac::I=1, 
    ) where {I <: Integer}
    n_atoms = length(atoms)
    if random == true
        ixs = rand(1:n_atoms, n_vac)
    end
    atoms_vac = [atoms[i] for i in 1:n_atoms if !(i in ixs)]
    coords_vac = [coords[i] for i in 1:n_atoms if !(i in ixs)]
    return atoms_vac, coords_vac
end

function make_supercell(atoms::Array, coords::Array, 
        box::Array, box_size::Diagonal; nx::I=1, ny::I=1,
        nz::I=1) where {I <: Integer}
    @assert (nx > 0) & (ny > 0) & (nz > 0) 
    sc_atoms = []
    sc_coords = []
    sc_box = box
    sc_box_size = box_size * Diagonal([nx, ny, nz])
    n_atoms = length(coords)
    for i in 0:nx-1, j in 0:ny-1, k in 0:nz-1
        push!(sc_atoms,atoms)
        scale = Diagonal([i,j,k])
        shift = sum(sc_box*scale, dims=1)
        push!(sc_coords,[coords[l]+shift for l in 1:n_atoms])
    end
    sc_atoms = vcat(sc_atoms...)
    sc_coords = vcat(sc_coords...)
    return sc_atoms, sc_coords, sc_box, sc_box_size
end

function plot_crystal(atoms::Array, coords::Array; 
        default_color::String="blue",
        element_color_map::Dict=Dict{String,String}(),
        default_size::T=50,
        element_size_map::Dict=Dict{String,Any}()
    ) where T <: Real
    
    elements = Set([atom.name for atom in atoms])
    for element in elements
        if !haskey(element_color_map, element)
            element_color_map[element] = default_color
        end
        if !haskey(element_size_map, element)
            element_size_map[element] = default_size
        end
    end
    colors = [element_color_map[element] for element in elements]
    sizes = [element_size_map[element] for element in elements]

    x = [v[1] for v in coords]
    y = [v[2] for v in coords]
    z = [v[3] for v in coords]
    return @gif for i in range(0, stop=2π, length=100)
        scatter(x, y, z, camera=(10*(1+cos(i)),5),
            markersize=sizes, legend=false, 
            color=colors, aspect_ratio=:equal,
            xlabel=L"x", ylabel=L"y", zlabel=L"z",
            title=string(length(atoms), " atoms of: ", join(elements, ","))
        )
    end
end

mutable struct MinimalSimulationConfig
    atoms::Array
    box_size::Float32
    coords::Array
    neighbours::Array{Tuple{Int64,Int64}}
end

struct NeighbourFinder 
    nb_matrix::BitArray{2} # defines which atom pairs we'll be happy to check at all
    n_steps::Int
    dist_cutoff::Float32
    rcut2::Float32
end

NeighbourFinder(nb_matrix, n_steps, dist_cutoff) = NeighbourFinder(nb_matrix, n_steps, dist_cutoff, dist_cutoff^2)

function vector1D(c1::Real, c2::Real, box_size::Real)
    if c1 < c2
        return (c2 - c1) < (c1 - c2 + box_size) ? (c2 - c1) : (c2 - c1 - box_size)
    else
        return (c1 - c2) < (c2 - c1 + box_size) ? (c2 - c1) : (c2 - c1 + box_size)
    end
end

vector(c1, c2, box_size::Real) = vector1D.(c1, c2, box_size)

function simple_find_neighbours(s::MinimalSimulationConfig,
        nf::NeighbourFinder, step_n::Int;
        parallel::Bool=false, 
        x_shifts=[0], y_shifts=[0], z_shifts=[0] # factors by which the box will be shifted along each box vector
    )
    
    !iszero(step_n % nf.n_steps) && return
    neighbours = empty(s.neighbours)
    for i in 1:length(s.coords)
        ci = s.coords[i]
        for j in 1:length(s.coords)
            if i==j 
                continue
            end
            
            r2 = sum(abs2, vector(ci, s.coords[j], s.box_size))
            if r2 <= nf.rcut2 && nf.nb_matrix[j,i]
                push!(neighbours, (i,j))
            end                
        end
    end
    return neighbours
end

function get_distance_df(atoms, box_size, coords; 
        initial_neighbours=Tuple{Int,Int}[],
        dist_cutoff=2
    )
    s = MinimalSimulationConfig(atoms, box_size, coords, initial_neighbours)
    n_atoms = length(atoms)
    nb_matrix = trues(n_atoms,n_atoms)
    n_steps = 1
    nf = NeighbourFinder(nb_matrix, n_steps, dist_cutoff)
    idxs = simple_find_neighbours(s, nf, 1)
    rs = [sqrt(sum(abs2, vector(s.coords[i], s.coords[j], s.box_size)))
        for (i,j) in idxs
    ]
    rs_df = sort(combine(groupby(DataFrame("distances"=>rs),[:distances]), 
 nrow=>:count), [:distances])
    return rs_df, rs
end

function plot_distance_hist(rs::Array, title::String, cutoff)
    histogram(rs, xlabel=L"r", ylabel="Frequency", 
        title=string(title, ": euclidan (periodic) distance distribution (rcut ",cutoff,")"),
        bins=200,
    )
end

end # module
