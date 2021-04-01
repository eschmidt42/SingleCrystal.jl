using Crystal
using Test

@testset "bcc->unit cell->super cell->vacancy" begin
    element = "Fe"
    mass = 42.
    el2atom_map = Dict(element => Crystal.Atom(name=element, mass=mass))
    a = 2.9
    unitcell = Crystal.make_bcc_unitcell([element for _ in 1:2], a, el2atom_map)
    # test that unitcell.box approx 2.9*I(3)
    # test that edge_lengths approx 2.9*ones(3)
    @test length(unitcell.coords) == length(unitcell.atoms)
    @test length(unitcell.coords) == 2

    supercell = Crystal.make_supercell(unitcell, nx=2, ny=3, nz=4)
    @test length(supercell.coords) == length(supercell.atoms)
    @test length(supercell.coords) == 2*2*3*4
    # test supercell box x = nx * unitcell box x
    # test supercell box y = ny * unitcell box y
    # test supercell box z = nz * unitcell box z

    vac = Crystal.add_vacancies(supercell, ixs=[1])
    @test length(vac.coords) == length(supercell.coords) - 1
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