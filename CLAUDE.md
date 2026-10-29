# CLAUDE.md

Istruzioni per Claude Code su questo repository. Si lavora e si scrive in **italiano**: codice,
commenti, documenti e messaggi di commit.

## Cos'è

Morfogenetica: cellule con una macchina a stati che fanno crescere forme da sole regole locali.
Si sta riscrivendo il motore in **Julia**, partendo da zero, su fondamenta diverse dal vecchio
motore JavaScript (congelato in `riferimento/`).

Prima di toccare il codice leggi, in quest'ordine:

1. `docs/specifica-motore-julia.md` — cosa deve fare il motore e perché
2. `docs/risultati-negativi.md` — cosa è già stato provato e **non** funziona
3. `docs/algebra-delle-primitive-morfogenetiche.md` — il quaderno di ricerca
4. `prove/attesi.md` — i numeri del motore vecchio (riferimento, non bersaglio)

I documenti di riferimento sono quelli in `docs/`: le vecchie copie su claude.ai sono archiviate e
non si aggiornano. Lo stato del lavoro è in `docs/piano.md`.

## Le tre regole di lavoro (non negoziabili)

1. **Il motore si congela quando è validato.** Dopo il tag di validazione, una modifica a
   `motore/` richiede un *test fisico* che fallisce. Una forma che si comporta male non è un motivo
   per cambiare il motore: si ridisegna la primitiva.
2. **Ogni grandezza che una cellula usa è uno dei quattro canali**: dentro di sé, dalle vicine,
   dalla matrice, misura assoluta. Se una regola ha bisogno di altro, la regola non è ammissibile:
   non si aggiunge una scorciatoia.
3. **Ogni comportamento deve essere chiedibile a una cellula vera.** Se è un dispositivo del
   simulatore (come lo `stampo` o l'`impalcatura`), non entra.

## Regole di metodo

- Una regola di fase riceve solo la struttura di percezione, mai lo stato globale. Una scansione
  globale per **leggere** (l'osservatore) è ammessa; una per **decidere** no.
- Ogni meccanismo ha un interruttore; la batteria riporta cosa succede spegnendolo.
- Casualità sempre con seme: ogni corsa dev'essere riproducibile.
- Unità vere (µm, min, nN, Pa). Ogni parametro viene dalla letteratura, con la fonte, oppure è
  dichiarato libero.
- I risultati negativi diventano test di regressione.
- Non sviluppare il motore JavaScript in `riferimento/`: serve solo come riferimento e per il banco.

## Struttura

```
src/Morfogenetica.jl  il modulo: include i file di motore/ e primitive/
motore/               geometria (vertici e facce), energia, topologia, chimica, passata
primitive/            una regola per file, pura, con il suo predicato di fine
programmi/            le catene, come tabelle di fasi
prove/                numeri attesi; test/ contiene la batteria eseguibile
test/runtests.jl      la batteria
scripts/prove.sh      lancia la batteria
riferimento/          motore JavaScript CONGELATO
confronto/            lo stesso esperimento in un motore maturo
```

## Le prove

Non ci sono test automatici su GitHub: la batteria si lancia a mano dalle sessioni, **prima di
ogni commit che tocca `motore/`, `primitive/` o `programmi/`**, e sempre quando l'utente lo chiede.

```
scripts/prove.sh
```

Se una prova fallisce, lo si dice chiaramente all'utente: non si commenta, non si salta, non si
allenta la tolleranza per farla passare.

## Installare Julia

Serve Julia ≥ 1.10. Nelle sessioni cloud lo installa lo script di setup dell'ambiente; a mano:

```
curl -fsSL https://install.julialang.org | sh -s -- --yes
export PATH="$HOME/.juliaup/bin:$PATH"
```

La rete deve permettere `*.julialang.org` (installazione e server dei pacchetti, compresi quelli
regionali come `us-west.pkg.julialang.org`) e `storage.julialang.net` — attenzione, **.net**: il
server dei pacchetti rimanda lì. Se `storage.julialang.net` è bloccato, si può scaricare il
registro da GitHub con `JULIA_PKG_SERVER="" scripts/prove.sh`.

## Stato

Fase 0, 1 e 2 fatte. La meccanica del modello laterale passa le prove F1–F5 (F5 contro tyssue) ed
è **congelata** (`motore-meccanica-v1`, segnato dal commit "CONGELATO: MOTORE-MECCANICA-V1"): da
qui `motore/01`–`05` si toccano solo con un test fisico che fallisce. La batteria lo fa
rispettare: `prove/congelati.txt` tiene l'impronta SHA-256 di ogni file congelato e la prova
"file congelati" fallisce se cambia. Una modifica ammessa aggiorna l'impronta nello stesso
commit del test che l'ha resa necessaria. Moduli nuovi si aggiungono come file nuovi, ognuno congelato quando è validato. Congelato
anche `06_integratore.jl` (`motore-integratore-v1`: Newton per gli equilibri, Eulero implicito
per la dinamica; per le forme si usa questo, non `rilassa!`). Funziona con qualunque tessuto che
fornisca energia, gradiente e controlli: l'interfaccia è descritta in cima al file. Congelato anche `07_guscio.jl` (`motore-guscio-v1`: il guscio 3D; i cambi di topologia
ricostruiscono il tessuto con `tessuto_da_poligoni`). Congelato anche `08_scambi.jl`
(`motore-scambi-v1`: scambi di vicine registrati, solo quando è corta tutta la parete). Da fare:
estrusione, scambio al bordo, chimica. Stato in
`docs/piano.md`, sezioni 6–9. Ogni file nuovo in `motore/` va aggiunto a mano in
`src/Morfogenetica.jl`: la cache di Julia non vede i file nuovi in una cartella.
