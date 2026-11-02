# Prove lente e prove veloci.
#
# `scripts/prove.sh` lancia tutta la batteria (circa 10 minuti); `scripts/prove.sh veloci` salta
# i gruppi marcati `@lenta` (quelli che portano un tessuto 3D all'equilibrio molte volte con
# Newton a matrice piena) e dice quanti ne ha saltati. Le veloci si lanciano durante il lavoro;
# la batteria intera prima di ogni commit che tocca `motore/` e prima di aprire o aggiornare una
# PR (CLAUDE.md).

const SOLO_VELOCI = get(ENV, "MORFOGENETICA_PROVE", "tutte") == "veloci"
const LENTE_SALTATE = Ref(0)

macro lenta(ex)
    return quote
        if SOLO_VELOCI
            LENTE_SALTATE[] += 1
            nothing
        else
            $(esc(ex))
        end
    end
end
