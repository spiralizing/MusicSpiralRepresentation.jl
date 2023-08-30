using MusicSpiralRepresentation
using Test

@testset "MusicSpiralRepresentation.jl" begin
    # Write your tests here.
    @test get_distance_to_keys(get_center_effect([60,64,67]))[1,1] == "C"
    #include test for reading csv/xml and 
end
