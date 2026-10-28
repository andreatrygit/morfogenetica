# CONGELATO

Questo è il motore JavaScript del lavoro precedente. **Non si sviluppa più.** Resta qui per due
ragioni: fare da riferimento numerico (vedi `../prove/attesi.md`) e far girare il banco
interattivo nel browser, che è utile per mostrare le cose.

Perché è stato abbandonato: le cellule sono sfere con una sola lunghezza di riposo per contatto,
quindi non possono diventare un cuneo e non piegano un foglio; le forze sono inventate invece di
discendere da un'energia; non c'è chimica. Il resto è in `../docs/specifica-motore-julia.md`,
sezione 1.

Contiene anche tre regole che **sembrano locali e non lo sono** (un contatore globale di quiete,
un controllo di stabilità su tutte le cellule, l'orologio del simulatore usato come orologio della
cellula). Sono elencate nella specifica, sezione 8: non vanno riportate.

## Cosa c'è

- `macchine.js` — libreria delle macchine a stati e costanti condivise
- `motore3d.js` — il motore, 2D e 3D (interruttore `piano`)
- `suite.js` — batteria delle quattro forme, con ferita
- `ablaz.js` — ablazione dei meccanismi, uno per uno
- `p2_gastrula.js` / `p3_gastrula.js` — la catena di gastrulazione in 2D e in 3D
- `p_simmetria.js` — la prova della legge di simmetria sulle elezioni
- `diag2.js` — profili radiali (frazione libera per fascia di raggio)
- `vedi.py` — vista ASCII di una sezione, da un dump JSON
- `primitive3d.js`, `simulatore.js`, `banco3d.html` — lo strato vecchio delle primitive (contiene
  `stampo`, `spingi`, `piatto`, `mangia`: tutti espulsi) e il banco nel browser

## Come si esegue

```
node -e "globalThis.window=globalThis; for(const f of ['macchine.js','motore3d.js','suite.js'])
  eval(require('fs').readFileSync(f,'utf8'));"
```

Gli script leggono i parametri dall'ambiente (`RR`, `RIN`, `KC`, `NP`, `MODI`, …): i nomi stanno
in testa a ciascun file.
