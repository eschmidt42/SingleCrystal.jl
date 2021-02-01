using Crystal
using Test

@testset "fcc" begin
    element = "Fe"
    a = 1
    atoms, coords, box, box_size, box_vectors = Crystal.make_fcc_unitcell(element, a=a)
    
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
    element = "Fe"
    a = 1
    atoms, coords, box, box_size, box_vectors = Crystal.make_bcc_unitcell(element, a=a)
    
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