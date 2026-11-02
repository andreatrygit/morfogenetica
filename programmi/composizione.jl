# Composizione delle fasi (specifica, sezione 5): `;` è la sequenza di fasi della tabella;
# `insieme` è ‖, due regole attive nella stessa fase.
#
# La legge quantitativa del conflitto non c'è ancora (specifica, sezione 5: "è esattamente il
# caso di cui ci manca la legge"). Quindi `insieme` accetta solo regole che non si
# contraddicono: se toccano la stessa intenzione (la stessa tensione, il volume bersaglio), la
# stessa voce dello stato o della dichiarazione con valori diversi, o chiedono cose diverse, si
# rifiuta con un errore invece di inventare una regola. Secrezioni e deposizioni si sommano:
# sono quantità.

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

function _unisci_intenzioni(riposo, a, b)
    campi = fieldnames(Intenzioni3D)
    valori = map(campi) do f
        r, x, y = getfield(riposo, f), getfield(a, f), getfield(b, f)
        (x != r && y != r && x != y) &&
            throw(ArgumentError("insieme: conflitto sull'intenzione $f ($x e $y), la legge manca"))
        x != r ? x : y
    end
    return Intenzioni3D(valori...)
end

function _unisci(p::Percezione, u::Uscita, v::Uscita)
    richiesta, piano = u.richiesta, u.piano
    if v.richiesta != :nessuna
        (richiesta != :nessuna && (richiesta != v.richiesta || !isequal(piano, v.piano))) &&
            throw(ArgumentError("insieme: richieste in conflitto ($richiesta e $(v.richiesta))"))
        richiesta, piano = v.richiesta, v.piano
    end
    return Uscita(_unisci_dizionari(u.stato, v.stato, "stato"),
                  _unisci_intenzioni(p.intenzioni, u.intenzioni, v.intenzioni),
                  _somma_dizionari(u.secrezioni, v.secrezioni),
                  _unisci_dizionari(u.dichiara, v.dichiara, "dichiara");
                  deposita = _somma_dizionari(u.deposita, v.deposita), richiesta, piano)
end

"""
    insieme(a, b; nome)

La fase con le regole delle fasi `a` e `b` attive insieme (‖). Finisce quando sono finite
tutte e due. Se le due regole si contraddicono, la regola della fase solleva un errore.
"""
function insieme(a::Fase, b::Fase; nome)
    regola = (p, par) -> _unisci(p, a.regola(p, par.a), b.regola(p, par.b))
    fine = (p, par) -> a.fine(p, par.a) && b.fine(p, par.b)
    return Fase(nome, regola, (a = a.parametri, b = b.parametri), fine)
end

"""
    con_fine(f, fine)

La stessa fase `f` con il predicato di fine `fine(p, parametri)`.
"""
con_fine(f::Fase, fine) = Fase(f.nome, f.regola, f.parametri, fine)
