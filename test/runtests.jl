# Batteria delle prove. Si lancia con `scripts/prove.sh` (vedi CLAUDE.md).
using Test
using Random
using Morfogenetica

include("previsioni.jl")

@testset "Morfogenetica" begin
    @testset "il pacchetto si carica" begin
        @test isdir(joinpath(Morfogenetica.RADICE, "docs"))
    end
    include("meccanica.jl")
    include("fisica_laterale.jl")
    include("contatto.jl")
    include("fisica_confronto.jl")
    include("integratore.jl")
    include("congelati.jl")
end
