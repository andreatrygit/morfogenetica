# Ogni primitiva deve essere plausibile (CLAUDE.md, regola 3; specifica, sezione 7): il suo file
# dice il meccanismo cellulare che la realizza e una fonte che lo documenta. Una primitiva senza
# meccanismo difendibile non entra.
#
# Commit rosso: la prova è marcata @test_broken finché ogni file non ha le due sezioni.

function prova_ogni_primitiva_ha_meccanismo_e_fonte()
    cartella = joinpath(Morfogenetica.RADICE, "primitive")
    file = filter(f -> endswith(f, ".jl"), readdir(cartella))
    isempty(file) && return false
    for f in file
        righe = readlines(joinpath(cartella, f))
        ha(prefisso) = any(r -> startswith(r, prefisso) && length(strip(r)) > length(prefisso) + 2, righe)
        (ha("# Meccanismo:") && ha("# Fonte:")) || return false
    end
    return true
end

@testset "le primitive sono plausibili" begin
    @test_broken prova_ogni_primitiva_ha_meccanismo_e_fonte()
end
