# Morfogenetica

Una disciplina inventata: **N sfere, ognuna con una macchina a stati, che fanno crescere forme
da sole regole locali.** Nessun seme, nessun piano centrale, nessuna cellula che sappia cosa
si sta costruendo.

Questo repository segna un **ripartire**. Il lavoro precedente (un motore in JavaScript con
cellule sferiche) ha prodotto forme che funzionano, un'algebra delle primitive, e una collezione
di risultati negativi che è la parte più preziosa. Ma il motore è diventato un sedimento di
correzioni fatte forma per forma, e va riscritto in Julia su fondamenta diverse.

## Da dove partire

1. **`docs/specifica-motore-julia.md`** — la specifica della riscrittura. Si legge per prima:
   dice cosa deve esserci nel motore nuovo, perché, e quale fallimento misurato motiva ogni
   decisione. Contiene anche come un programma si compila nelle primitive di una cellula e la
   descrizione di una passata del motore.
2. **`docs/algebra-delle-primitive-morfogenetiche.md`** — il quaderno: le primitive, le regole
   di buona formazione, gli enunciati, i lavori vicini verificati, e i risultati negativi con i
   loro numeri. È il documento di ricerca; la specifica ne è la conseguenza ingegneristica.
3. **`docs/risultati-negativi.md`** — l'estratto operativo: ciò che è stato provato e **non**
   funziona, da far valere come test di regressione nel motore nuovo.
4. **`prove/attesi.md`** — i numeri che il motore vecchio produce, da usare come riferimento
   (non come bersaglio: le unità cambieranno).

## Le tre regole di lavoro

Vengono dall'errore di metodo che ha reso necessaria la riscrittura.

1. **Il motore si congela quando è validato.** Dopo, una modifica richiede un *test fisico* che
   falla — non una forma che si comporta male.
2. **Ogni grandezza che una cellula usa va dichiarata** come uno dei quattro canali (dentro di
   sé, dalle vicine, dalla matrice, misura assoluta). Se una regola ha bisogno di qualcosa che
   non c'è, non si aggiunge una scorciatoia: si dichiara che la regola non è ammissibile.
3. **Ogni comportamento deve essere chiedibile a una cellula vera.** Se la risposta alla domanda
   "cosa chiederesti a una cellula o a un organoide per ottenerlo?" è "nulla, è un dispositivo
   del simulatore", il comportamento non entra.

## Struttura

```
docs/         specifica, algebra, risultati negativi
motore/       Julia: geometria (vertici e facce), energia, topologia, chimica, passata
primitive/    Julia: una regola per file, pura, con il suo predicato di fine
programmi/    le catene, come tabelle di fasi
prove/        batteria con i numeri attesi; i risultati negativi come regressioni
riferimento/  il motore JavaScript CONGELATO + il banco nel browser
confronto/    lo stesso esperimento in un motore maturo, per validare la fisica
```

## Primo traguardo

Non una forma: un **test fisico**. Un foglio epiteliale piatto, una fila di cellule che alza la
tensione dell'anello apicale, e la curvatura che ne risulta confrontata con la previsione
analitica del modello a vertici. Se quello torna, la fisica tiene e si può ricominciare a
costruire forme — questa volta senza toccare il motore.

## I documenti vivi

I due documenti principali sono anche documenti vivi su claude.ai, dove vengono aggiornati:

- algebra: https://claude.ai/artifact/28ddvGrMSXCBoKfiaG4bVi
- specifica: https://claude.ai/artifact/8y9AazU8HEowgy3oKMV482

Le copie in `docs/` sono istantanee del 27 ottobre 2026. Se divergono, il documento vivo vince.
