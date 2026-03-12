using MusicSpiralRepresentation
using Test

@testset "MusicSpiralRepresentation.jl" begin
    @test get_distance_to_keys(get_center_effect([60,64,67]))[1,1] == "C"
end
