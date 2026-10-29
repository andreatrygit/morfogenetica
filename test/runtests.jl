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
    include("guscio.jl")
    include("scambi.jl")
    include("estrusione.jl")
    include("scutoidi.jl")
    include("fisica_confronto_3d.jl")
    include("congelati.jl")
end
