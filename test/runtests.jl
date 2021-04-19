using Crystal
using Test
using LinearAlgebra
using DataFrames
using JSON

@testset "fcc/bcc unit cell->super cell->vacancy" begin
    element = "Fe"
    mass = 42.
    el2atom_map = Dict(element => Crystal.Atom(name=element, mass=mass))
    a = 2.9
    nx, ny, nz = 2, 3, 4
    for b in ["fcc", "bcc"]
        if b == "fcc"
            unitcell = Crystal.make_fcc_unitcell(element, a)
            num = 4
        elseif b == "bcc"
            unitcell = Crystal.make_bcc_unitcell(element, a)
            num = 2
        end
        @testset "$(b)" begin
            
            @testset "unitcell: box" begin
                @test unitcell.box.M ≈ a*I(3)
                @test unitcell.edge_lengths ≈ a*ones(3)
            end
            @testset "unitcell: positions & atoms" begin
                @test length(unitcell.positions) == length(unitcell.atoms)
                @test length(unitcell.positions) == num
            end

            supercell = Crystal.make_supercell(unitcell, nx=nx, ny=ny, nz=nz)
            @testset "supercell: box" begin
                @test supercell.box.M[1,1] ≈ nx * unitcell.box.M[1,1]
                @test supercell.box.M[2,2] ≈ ny * unitcell.box.M[2,2]
                @test supercell.box.M[3,3] ≈ nz * unitcell.box.M[3,3]
            end
            @testset "supercell: positions & atoms" begin
                @test length(supercell.positions) == length(supercell.atoms)
                @test length(supercell.positions) == num*2*3*4
            end
            
            @testset "vacancy" begin
                vac = Crystal.add_vacancies(supercell, ixs=[1])
                @test length(vac.positions) == length(supercell.positions) - 1
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
    for i in 1:length(cell.positions)
        ci = cell.positions[i]
        for j in 1:length(cell.positions)
            if i==j 
                continue
            end
            
            r2 = sum(abs2, Crystal.vector(ci, cell.positions[j], cell.edge_lengths))
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
        unitcell = Crystal.make_bcc_unitcell(el, 1.,)
        supercell = Crystal.make_supercell(unitcell, nx=3, ny=3, nz=3)
        rs_df, rs = get_distance_df(supercell, dist_cutoff=d)
        @test rs_df.distances[1] ≈ sqrt((sqrt(2)/2)^2 + 1/2^2)
        @test rs_df.distances[2] ≈ 1
        @test rs_df.distances[3] ≈ sqrt(2)
        @test rs_df.distances[4] ≈ sqrt((sqrt(2)/2)^2 + (3/2)^2)
        @test rs_df.distances[5] ≈ sqrt(sqrt(2)^2 + 1^2)
    end

    @testset "fcc" begin
        unitcell = Crystal.make_fcc_unitcell(el, 1.)
        supercell = Crystal.make_supercell(unitcell, nx=3, ny=3, nz=3)
        rs_df, rs = get_distance_df(supercell, dist_cutoff=d)
        @test rs_df.distances[1] ≈ sqrt(1^2+1^2)/2
        @test rs_df.distances[2] ≈ 1
        @test rs_df.distances[3] ≈ sqrt(1^2+(sqrt(2)/2)^2)
        @test rs_df.distances[4] ≈ sqrt(1^2+1^2)
        @test rs_df.distances[5] ≈ sqrt(3^2+1^2)/2
    end
end

crystal_specs = Dict(
    # NaCl structure
    "NaCl" => Dict(
        "symbols" => ["Na", "Cl"],
        "basis" => [[0. 0. 0.], [.5 .5 .5]], # scaled coordinates
        "nr" => 225,
        "setting" => 1,
        "a_direction" => [1.; 0.; 0.],
        "ab_normal" => [0.; 0.; 1.],
        "cellpar" => [5.64, 5.64, 5.64, 90, 90, 90]
    ),
    # Al fcc structure
    "Al_fcc" => Dict(
        "symbols" => ["Al"],
        "basis" => [[0. 0. 0.],], # scaled coordinates
        "nr" => 225,
        "setting" => 1,
        "a_direction" => [1.; 0.; 0.],
        "ab_normal" => [0.; 0.; 1.],
        "cellpar" => [4.05, 4.05, 4.05, 90, 90, 90]
    ),
    # Fe bcc structure
    "Fe_bcc" => Dict(
        "symbols" => ["Fe"],
        "basis" => [[0. 0. 0.],], # scaled coordinates
        "nr" => 229,
        "setting" => 1,
        "a_direction" => [1.; 0.; 0.],
        "ab_normal" => [0.; 0.; 1.],
        "cellpar" => [2.87, 2.87, 2.87, 90, 90, 90]
    ),
    # Mg hcp structure
    "Mg_hcp" => Dict(
        "symbols" => ["Mg"],
        "basis" => [[1/3 2/3 3/4],], # scaled coordinates
        "nr" => 194,
        "setting" => 1,
        "a_direction" => [1.; 0.; 0.],
        "ab_normal" => [0.; 0.; 1.],
        "cellpar" => [3.21, 3.21, 5.21, 90, 90, 120]
    ),
    # Diamond structure
    "Diamond" => Dict(
        "symbols" => ["C"],
        "basis" => [[0. 0. 0.],], # scaled coordinates
        "nr" => 227,
        "setting" => 1,
        "a_direction" => [1.; 0.; 0.],
        "ab_normal" => [0.; 0.; 1.],
        "cellpar" => [3.57, 3.57, 3.57, 90, 90, 90]
    ),
    # Rutile structure
    "Rutile" => Dict(
        "symbols" => ["Ti", "O"],
        "basis" => [[0. 0. 0.], [.3 .3 0.]], # scaled coordinates
        "nr" => 136,
        "setting" => 1,
        "a_direction" => [1.; 0.; 0.],
        "ab_normal" => [0.; 0.; 1.],
        "cellpar" => [4.6, 4.6, 2.95, 90, 90, 90]
    ),
    # CoSb3 skudderudite
    "Skudderudite" => Dict(
        "symbols" => ["Co", "Sb"],
        "basis" => [[.25 .25 .25], [0. .335 .158]], # scaled coordinates
        "nr" => 204,
        "setting" => 1,
        "a_direction" => [1.; 0.; 0.],
        "ab_normal" => [0.; 0.; 1.],
        "cellpar" => [9.04, 9.04, 9.04, 90, 90, 90]
    )
)

function positions_match_ase(crystal::Crystal.Cell, ase_crystal::Dict)
    return all([p0 ≈ p1 for (p0,p1) in zip(ase_crystal["positions"],crystal.positions)])
end

function cell_matches_ase(crystal::Crystal.Cell, ase_crystal::Dict)
    return crystal.box.M ≈ ase_crystal["cell"]
end

@testset "various single crystals" begin

    json_ase_crystals = Crystal.load_refs()

    for name in keys(crystal_specs)
        @testset "$(name)" begin
            nr = crystal_specs[name]["nr"]
            setting = crystal_specs[name]["setting"]
            basis = crystal_specs[name]["basis"]
            symbols = crystal_specs[name]["symbols"]
            a_direction = crystal_specs[name]["a_direction"]
            ab_normal = crystal_specs[name]["ab_normal"]
            cellpar = crystal_specs[name]["cellpar"]
    
            crystal = Crystal.make_unitcell(basis, symbols, nr, setting, cellpar,
                                   a_direction=a_direction, ab_normal=ab_normal)
            
            ase_crystal = Crystal.parse_json_crystal(json_ase_crystals[name])
            
            @testset "positions match" begin 
                @test positions_match_ase(crystal, ase_crystal)
            end
            @testset "cell match" begin
                @test cell_matches_ase(crystal, ase_crystal)
            end
        end
    end
end