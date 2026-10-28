# Batteria delle prove. Si lancia con `scripts/prove.sh` (vedi CLAUDE.md).
using Test
using Morfogenetica

@testset "Morfogenetica" begin
    @testset "il pacchetto si carica" begin
        @test isdefined(Morfogenetica, :RADICE)
        @test isdir(joinpath(Morfogenetica.RADICE, "docs"))
    end
end
