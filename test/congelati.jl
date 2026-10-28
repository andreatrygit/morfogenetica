# Il congelamento del motore fatto rispettare dalla batteria (prove/congelati.txt).

using SHA

@testset "file congelati" begin
    elenco = joinpath(Morfogenetica.RADICE, "prove", "congelati.txt")
    voci = [split(r) for r in readlines(elenco) if !isempty(strip(r)) && !startswith(r, "#")]
    @test !isempty(voci)
    for (impronta, file) in voci
        percorso = joinpath(Morfogenetica.RADICE, file)
        @test isfile(percorso)
        attuale = bytes2hex(open(sha256, percorso))
        attuale == impronta || @warn "file congelato modificato: $file (vedi prove/congelati.txt)"
        @test attuale == impronta
    end
end
