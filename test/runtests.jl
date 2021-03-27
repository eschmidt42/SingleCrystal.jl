using Crystal
using Test

masses = Dict("V" => 50.9415, "Nb" => 92.9064, "Ta" => 180.9479,
              "Cr" => 51.996, "Mo" => 95.94, "W" => 183.85,
              "Fe" => 55.847)

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