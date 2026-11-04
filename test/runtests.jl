# Batteria delle prove. Si lancia con `scripts/prove.sh` (tutte) o `scripts/prove.sh veloci`
# (salta i gruppi marcati @lenta; vedi test/lente.jl e CLAUDE.md).
using Test
using Random
using Morfogenetica

include("previsioni.jl")
include("lente.jl")

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
    include("operazioni.jl")
    include("operazioni_periodiche.jl")
    include("chimica.jl")
    include("programma.jl")
    include("popolazione.jl")
    include("invaginazione.jl")
    include("composizione.jl")
    include("primitive.jl")
    include("organizzatore.jl")
    include("induci.jl")
    include("prestazioni.jl")
    include("estensioni.jl")
    include("sparso.jl")
    include("gastrula.jl")
    include("fisica_confronto_3d.jl")
    include("congelati.jl")
end

SOLO_VELOCI && @info "Prove veloci: saltati $(LENTE_SALTATE[]) gruppi lenti. La batteria intera è scripts/prove.sh."
