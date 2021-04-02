using Crystal
using Test
using LinearAlgebra
using DataFrames

@testset "unit cell->super cell->vacancy" begin
    element = "Fe"
    mass = 42.
    el2atom_map = Dict(element => Crystal.Atom(name=element, mass=mass))
    a = 2.9
    nx, ny, nz = 2, 3, 4
    for b in ["fcc", "bcc"]
        if b == "fcc"
            elements = [element for _ in 1:4]
        elseif b == "bcc"
            elements = [element for _ in 1:2]
        end
        @testset "$(b)" begin
            
            unitcell = Crystal.make_unitcell(elements, a, el2atom_map, b)
            @testset "unitcell: box" begin
                @test unitcell.box.M ≈ a*I(3)
                @test unitcell.edge_lengths ≈ a*ones(3)
            end
            @testset "unitcell: coords & atoms" begin
                @test length(unitcell.coords) == length(unitcell.atoms)
                @test length(unitcell.coords) == length(elements)
            end

            supercell = Crystal.make_supercell(unitcell, nx=nx, ny=ny, nz=nz)
            @testset "supercell: box" begin
                @test supercell.box.M[1,1] ≈ nx * unitcell.box.M[1,1]
                @test supercell.box.M[2,2] ≈ ny * unitcell.box.M[2,2]
                @test supercell.box.M[3,3] ≈ nz * unitcell.box.M[3,3]
            end
            @testset "supercell: coords & atoms" begin
                @test length(supercell.coords) == length(supercell.atoms)
                @test length(supercell.coords) == length(elements)*2*3*4
            end
            
            @testset "vacancy" begin
                vac = Crystal.add_vacancies(supercell, ixs=[1])
                @test length(vac.coords) == length(supercell.coords) - 1
            end
        end
    end
end

struct NeighbourFinder{T} 
    nb_matrix::BitArray{2} # defines which atom pairs we'll be happy to check at all
    dist_cutoff::T
    rcut2::T
end

NeighbourFinder(nb_matrix, dist_cutoff) = NeighbourFinder{typeof(dist_cutoff)}(nb_matrix, dist_cutoff, dist_cutoff^2)

function find_neighbours(cell::Crystal.Cell, nf::NeighbourFinder)
    neighbours = []
    rs = []
    for i in 1:length(cell.coords)
        ci = cell.coords[i]
        for j in 1:length(cell.coords)
            if i==j 
                continue
            end
            
            r2 = sum(abs2, Crystal.vector(ci, cell.coords[j], cell.edge_lengths))
            if r2 <= nf.rcut2 && nf.nb_matrix[j,i]
                push!(neighbours, (i,j))
                push!(rs, sqrt(r2))
            end                
        end
    end
    return neighbours, rs
end

function get_distance_df(cell::Crystal.Cell; dist_cutoff::Real=2)
    n_atoms = length(cell.atoms)
    nb_matrix = trues(n_atoms,n_atoms)
    nf = NeighbourFinder(nb_matrix, dist_cutoff)
    idxs, rs = find_neighbours(cell, nf)
    rs_df = sort(combine(groupby(DataFrame("distances"=>rs),[:distances]), 
 nrow=>:count), [:distances])
    return rs_df, rs
end

@testset "Neighbours" begin
    
    d = 2.
    el = "Fe"
    el2atom_map = Dict(el => Crystal.Atom(name=el, mass=1))

    @testset "bcc" begin
        unitcell = Crystal.make_bcc_unitcell([el for _ in 1:2], 1., el2atom_map)
        supercell = Crystal.make_supercell(unitcell, nx=3, ny=3, nz=3)
        rs_df, rs = get_distance_df(supercell, dist_cutoff=d)
        @test rs_df.distances[1] ≈ sqrt((sqrt(2)/2)^2 + 1/2^2)
        @test rs_df.distances[2] ≈ 1
        @test rs_df.distances[3] ≈ sqrt(2)
        @test rs_df.distances[4] ≈ sqrt((sqrt(2)/2)^2 + (3/2)^2)
        @test rs_df.distances[5] ≈ sqrt(sqrt(2)^2 + 1^2)
    end

    @testset "fcc" begin
        unitcell = Crystal.make_fcc_unitcell([el for _ in 1:4], 1., el2atom_map)
        supercell = Crystal.make_supercell(unitcell, nx=3, ny=3, nz=3)
        rs_df, rs = get_distance_df(supercell, dist_cutoff=d)
        @test rs_df.distances[1] ≈ sqrt(1^2+1^2)/2
        @test rs_df.distances[2] ≈ 1
        @test rs_df.distances[3] ≈ sqrt(1^2+(sqrt(2)/2)^2)
        @test rs_df.distances[4] ≈ sqrt(1^2+1^2)
        @test rs_df.distances[5] ≈ sqrt(3^2+1^2)/2
    end
end


# sanity check first, second and higher nearest neighbours for edge lengths = 1

#=
masses = Dict("V" => 50.9415, "Nb" => 92.9064, "Ta" => 180.9479,
              "Cr" => 51.996, "Mo" => 95.94, "W" => 183.85,
              "Fe" => 55.847)

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

@testset "vacancy" begin
    el = "Fe"
    a = 1
    el2atom_map = Dict(el => Crystal.Atom(name=el, mass=masses[el]))
    elements = [el for _ in 1:4]
    atoms, coords, box, box_size, box_vectors = Crystal.make_fcc_unitcell(elements, a=a, el2atom_map=el2atom_map)
    vac_atoms, vac_coords = Crystal.add_vacancies(atoms, coords, ixs=[1])
    @test length(vac_atoms) == length(vac_coords)
    @test length(vac_atoms) == length(atoms) - 1
end

@testset "fcc" begin
    el = "Fe"
    a = 1
    el2atom_map = Dict(el => Crystal.Atom(name=el, mass=masses[el]))
    elements = [el for _ in 1:4]
    atoms, coords, box, box_size, box_vectors = Crystal.make_fcc_unitcell(elements, a=a, el2atom_map=el2atom_map)
    
    sc_atoms, sc_coords, sc_box, sc_box_size = Crystal.make_supercell(atoms, coords, box, box_size, nx=3, ny=3,
        nz=3);
    @test length(sc_atoms) == length(sc_coords)
    
    dist_cutoff = 2
    rs_df, rs = Crystal.get_distance_df(sc_atoms, sc_box_size[1,1], sc_coords, dist_cutoff=dist_cutoff)
    @test rs_df.distances[1] ≈ sqrt(1^2+1^2)/2
    @test rs_df.distances[2] ≈ 1
    @test rs_df.distances[3] ≈ sqrt(1^2+(sqrt(2)/2)^2)
    @test rs_df.distances[4] ≈ sqrt(1^2+1^2)
    @test rs_df.distances[5] ≈ sqrt(3^2+1^2)/2
end

@testset "bcc" begin
    el = "Fe"
    a = 1
    el2atom_map = Dict(el => Crystal.Atom(name=el, mass=masses[el]))
    elements = [el for _ in 1:2]
    atoms, coords, box, box_size, box_vectors = Crystal.make_bcc_unitcell(elements, a=a, el2atom_map=el2atom_map)
    
    sc_atoms, sc_coords, sc_box, sc_box_size = Crystal.make_supercell(atoms, coords, box, box_size, nx=3, ny=3,
        nz=3);
    @test length(sc_atoms) == length(sc_coords)
    
    dist_cutoff = 2
    rs_df, rs = Crystal.get_distance_df(sc_atoms, sc_box_size[1,1], sc_coords, dist_cutoff=dist_cutoff)
    @test rs_df.distances[1] ≈ sqrt((sqrt(2)/2)^2 + 1/2^2)
    @test rs_df.distances[2] ≈ 1
    @test rs_df.distances[3] ≈ sqrt(2)
    @test rs_df.distances[4] ≈ sqrt((sqrt(2)/2)^2 + (3/2)^2)
    @test rs_df.distances[5] ≈ sqrt(sqrt(2)^2 + 1^2)
end
=#