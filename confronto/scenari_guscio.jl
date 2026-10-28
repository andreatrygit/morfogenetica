# Scenari del confronto 3D con tyssue (prova F7): definiti una volta sola, qui.
# Usato da confronto/esporta_guscio.jl e da test/fisica_confronto_3d.jl.

const P_GUSCIO_CONFRONTO = ParametriMeccanici(K = 1e-3, η = 1.0, lunghezza_minima = 0.05,
    k_contatto = 0.0, d_contatto = 0.5)

# Pezzo di foglio con 3 anelli attorno alla cellula centrale (37 cellule), bordo bloccato.
function scenario_guscio(nome)
    R = 3
    n = 1 + 3R * (R + 1)
    if nome == "anello_centrale"
        # tensioni uguali, la cellula centrale stringe l'anello apicale
        ints = fill(Intenzioni3D(1.0, 1.0, 1.0, 0.0, 1000.0), n)
        ints[1] = Intenzioni3D(1.0, 1.0, 1.0, 2.0, 1000.0)
    elseif nome == "placca"
        # apicale più teso del basale ovunque; la cellula centrale e il primo anello stringono
        ints = fill(Intenzioni3D(1.2, 0.8, 1.0, 0.0, 1000.0), n)
        for c in 1:7
            ints[c] = Intenzioni3D(1.2, 0.8, 1.0, 1.5, 1000.0)
        end
    else
        throw(ArgumentError("scenario sconosciuto: $nome"))
    end
    return foglio_esagonale(ints, R; lato = 5.0, altezza = 15.0, bordo_fisso = true)
end

const SCENARI_GUSCIO = ("anello_centrale", "placca")
