# Composizione delle fasi (specifica, sezione 5): `;` è la sequenza di fasi della tabella;
# `insieme` è ‖, due regole attive nella stessa fase.
#
# La legge quantitativa del conflitto non c'è ancora (specifica, sezione 5: "è esattamente il
# caso di cui ci manca la legge"). Quindi `insieme` accetta solo fasi che non si contraddicono,
# e lo decide quando compone, dalle dichiarazioni (`Fase(...; scrive, stato, dichiara)`): due
# fasi che scrivono la stessa intenzione, la richiesta, la stessa voce dello stato interno o
# della dichiarazione alle vicine si rifiutano, e così le fasi non dichiarate, invece di
# inventare una regola. Secrezioni e deposizioni si sommano: sono quantità.
#
# Perché alla composizione e non alla corsa. Un controllo alla corsa (un registro delle
# scritture, versioni, lock) vede solo i rami eseguiti: una regola che scrive Λa solo sopra
# soglia rivela il conflitto solo quando qualche cellula supera la soglia, magari mai nelle
# prove. Le dichiarazioni lo vedono per ogni cellula. Il registro esatto delle scritture serve
# lo stesso, alla passata, per verificare che le dichiarazioni non mentano: è l'uscita a delta
# (13_programma.jl). Un lock deciderebbe un ordine, cioè "vince l'ultimo": proprio la legge
# che manca. Fra cellule diverse non ci sono conflitti per costruzione: ognuna scrive solo il
# proprio stato, e tutte leggono la passata precedente.

export insieme, con_fine

# Le uscite sono delta su voci disgiunte (le dichiarazioni lo garantiscono): si uniscono; le
# quantità si sommano.
function _disgiunti(a, b, cosa)
    comuni = intersect(keys(a), keys(b))
    isempty(comuni) || throw(ArgumentError("insieme: tutte e due scrivono $cosa $(collect(comuni))"))
    return merge(a, b)
end

function _unisci(u::Uscita, v::Uscita)
    (u.richiesta !== nothing && v.richiesta !== nothing) &&
        throw(ArgumentError("insieme: tutte e due scrivono la richiesta"))
    richiesta, piano = v.richiesta === nothing ? (u.richiesta, u.piano) : (v.richiesta, v.piano)
    return Uscita(; intenzioni = _disgiunti(u.intenzioni, v.intenzioni, "le intenzioni"),
                  stato = _disgiunti(u.stato, v.stato, "lo stato"),
                  dichiara = _disgiunti(u.dichiara, v.dichiara, "la dichiarazione"),
                  secrezioni = mergewith(+, u.secrezioni, v.secrezioni),
                  deposita = mergewith(+, u.deposita, v.deposita), richiesta, piano)
end

"""
    insieme(a, b; nome)

La fase con le regole delle fasi `a` e `b` attive insieme (‖). Finisce quando sono finite
tutte e due e scrive ciò che scrivono le due. Rifiuta (ArgumentError) le fasi non dichiarate e
due fasi che scrivono la stessa intenzione, la richiesta, la stessa voce di stato o di
dichiarazione.
"""
function insieme(a::Fase, b::Fase; nome)
    for f in (a, b)
        f.scrive === nothing &&
            throw(ArgumentError("insieme: la fase $(f.nome) non dichiara che cosa scrive"))
    end
    for (cosa, x, y) in (("", setdiff(a.scrive, (:secrezioni, :deposita)), b.scrive),
                         ("stato ", a.scrive_stato, b.scrive_stato),
                         ("dichiara ", a.scrive_dichiara, b.scrive_dichiara))
        comuni = intersect(x, y)
        isempty(comuni) ||
            throw(ArgumentError("insieme: $(a.nome) e $(b.nome) scrivono entrambe $cosa" *
                                "$(Tuple(comuni)); la legge del conflitto manca"))
    end
    regola = (p, par) -> _unisci(a.regola(p, par.a), b.regola(p, par.b))
    fine = (p, par) -> a.fine(p, par.a) && b.fine(p, par.b)
    return Fase(nome, regola, (a = a.parametri, b = b.parametri), fine;
                scrive = Tuple(union(a.scrive, b.scrive)),
                stato = Tuple(union(a.scrive_stato, b.scrive_stato)),
                dichiara = Tuple(union(a.scrive_dichiara, b.scrive_dichiara)))
end

"""
    con_fine(f, fine)

La stessa fase `f` con il predicato di fine `fine(p, parametri)`.
"""
con_fine(f::Fase, fine) = Fase(f.nome, f.regola, f.parametri, fine; scrive = f.scrive,
                               stato = something(f.scrive_stato, ()),
                               dichiara = something(f.scrive_dichiara, ()))
