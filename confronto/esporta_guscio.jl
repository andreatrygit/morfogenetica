# Esporta la geometria di partenza degli scenari 3D per tyssue (confronto/tyssue_guscio.py).
# Le tensioni sono già tradotte nelle convenzioni di tyssue:
# - tyssue dà a ogni cellula una copia propria di ogni faccia: una faccia laterale condivisa
#   ha due copie, ciascuna con metà della tensione media; una faccia di una cellula sola ha
#   la tensione di quella cellula;
# - LineTension vale λ·L/2 per mezzo-spigolo: per la nostra quota Λ/condivisione serve
#   λ = 2Λ/condivisione.
#
# Si lancia dalla radice:  julia --project=. confronto/esporta_guscio.jl

using Morfogenetica
include(joinpath(@__DIR__, "scenari_guscio.jl"))

const CARTELLA = joinpath(@__DIR__, "risultati")

for nome in SCENARI_GUSCIO
    t = scenario_guscio(nome)
    open(joinpath(CARTELLA, "guscio_$(nome)_vertici.csv"), "w") do io
        println(io, "vertice,x,y,z,fisso")
        for j in 1:length(t.x)÷3
            println(io, j, ",", t.x[3j-2], ",", t.x[3j-1], ",", t.x[3j], ",", Int(t.fissi[j]))
        end
    end
    open(joinpath(CARTELLA, "guscio_$(nome)_facce.csv"), "w") do io
        # una riga per copia di faccia: cellula, etichetta, tensione di superficie, vertici in
        # senso antiorario visti da fuori della cellula, tensioni di linea dei suoi spigoli
        println(io, "cellula,etichetta,tensione,vertici,tensioni_linea")
        for (c, cel) in enumerate(t.cellule), (fi, s) in zip(cel.facce, cel.segni)
            f = t.facce[fi]
            σ = Morfogenetica._tensione_faccia(t, f) / length(f.cellule)
            v = s == 1 ? f.vertici : reverse(f.vertici)
            λ = f.etichetta == apicale ?
                [2 * t.intenzioni[c].Λa / f.condivisione[k] for k in eachindex(v)] :
                zeros(length(v))
            println(io, c, ",", f.etichetta, ",", σ, ",", join(v, " "), ",", join(λ, " "))
        end
    end
    open(joinpath(CARTELLA, "guscio_$(nome)_cellule.csv"), "w") do io
        println(io, "cellula,V0,K")
        for c in eachindex(t.cellule)
            println(io, c, ",", t.intenzioni[c].V0, ",", P_GUSCIO_CONFRONTO.K)
        end
    end
    println(nome, ": ", length(t.cellule), " cellule, ", length(t.x) ÷ 3, " vertici, ",
            count(t.fissi), " bloccati")
end
