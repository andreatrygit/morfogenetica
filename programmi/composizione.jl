# Composizione delle fasi (specifica, sezione 5): `;` è la sequenza di fasi della tabella;
# `insieme` è ‖, due regole attive nella stessa fase.
#
# La legge quantitativa del conflitto non c'è ancora (specifica, sezione 5: "è esattamente il
# caso di cui ci manca la legge"). Quindi `insieme` accetta solo fasi che non si contraddicono,
# e lo decide quando compone, dalle dichiarazioni (`Fase(...; scrive)`): due fasi che scrivono
# la stessa intenzione (la stessa tensione, il volume bersaglio) o la richiesta si rifiutano, e
# così le fasi non dichiarate, invece di inventare una regola. Secrezioni e deposizioni si
# sommano: sono quantità. Lo stato interno e le dichiarazioni alle vicine si uniscono voce per
# voce, e due valori diversi per la stessa voce sono un errore alla corsa.

export insieme, con_fine

function _unisci_dizionari(a, b, cosa)
    u = copy(a)
    for (k, v) in b
        haskey(u, k) && u[k] != v &&
            throw(ArgumentError("insieme: conflitto su $cosa[$k] ($(u[k]) e $v), la legge manca"))
        u[k] = v
    end
    return u
end

_somma_dizionari(a, b) = mergewith(+, a, b)

# Ogni intenzione viene dalla fase che la dichiara (al più una), altrimenti resta a riposo.
function _unisci_intenzioni(riposo, a, sa, b, sb)
    valori = map(fieldnames(Intenzioni3D)) do f
        f in sa ? getfield(a, f) : f in sb ? getfield(b, f) : getfield(riposo, f)
    end
    return Intenzioni3D(valori...)
end

function _unisci(p::Percezione, u::Uscita, sa, v::Uscita, sb)
    richiesta, piano = :richiesta in sb ? (v.richiesta, v.piano) : (u.richiesta, u.piano)
    return Uscita(_unisci_dizionari(u.stato, v.stato, "stato"),
                  _unisci_intenzioni(p.intenzioni, u.intenzioni, sa, v.intenzioni, sb),
                  _somma_dizionari(u.secrezioni, v.secrezioni),
                  _unisci_dizionari(u.dichiara, v.dichiara, "dichiara");
                  deposita = _somma_dizionari(u.deposita, v.deposita), richiesta, piano)
end

"""
    insieme(a, b; nome)

La fase con le regole delle fasi `a` e `b` attive insieme (‖). Finisce quando sono finite
tutte e due e scrive ciò che scrivono le due. Rifiuta (ArgumentError) le fasi non dichiarate e
due fasi che scrivono la stessa intenzione o la richiesta.
"""
function insieme(a::Fase, b::Fase; nome)
    for f in (a, b)
        f.scrive === nothing &&
            throw(ArgumentError("insieme: la fase $(f.nome) non dichiara che cosa scrive"))
    end
    comuni = setdiff(intersect(a.scrive, b.scrive), (:secrezioni, :deposita))
    isempty(comuni) ||
        throw(ArgumentError("insieme: $(a.nome) e $(b.nome) scrivono entrambe $(Tuple(comuni)); " *
                            "la legge del conflitto manca"))
    sa, sb = a.scrive, b.scrive
    regola = (p, par) -> _unisci(p, a.regola(p, par.a), sa, b.regola(p, par.b), sb)
    fine = (p, par) -> a.fine(p, par.a) && b.fine(p, par.b)
    return Fase(nome, regola, (a = a.parametri, b = b.parametri), fine;
                scrive = Tuple(union(sa, sb)))
end

"""
    con_fine(f, fine)

La stessa fase `f` con il predicato di fine `fine(p, parametri)`.
"""
con_fine(f::Fase, fine) = Fase(f.nome, f.regola, f.parametri, fine; scrive = f.scrive)
