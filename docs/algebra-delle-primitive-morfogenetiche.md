# Algebra delle primitive morfogenetiche

Oct 23, 2026 · @Andrea

Una forma biologica si può scrivere come un programma: poche primitive locali, eseguite da ogni cellula sul proprio vicinato, combinate da cinque operatori. Questo documento definisce il modello, la semantica delle primitive e degli operatori, le leggi che dovrebbero valere e le regole che le prove nel Banco Morfogenetico hanno reso necessarie.

## 1. Perché un'algebra

Finora ogni forma è stata costruita a mano: un disco, un anello, un'ellisse, una sfera, un guscio, ognuno con le sue regole aggiunte quando una prova falliva. Due osservazioni suggeriscono che sotto c'è una struttura comune.

- **La macchina a stati non cambia.** Da quando il seme è sparito, tutte le forme usano la stessa logica: cresci se la figlia cade nel bersaglio, fermati, muori se sei fuori. Cambiano solo le coordinate che la cellula conosce di sé e il bersaglio.
- **Le stesse soluzioni ritornano.** Biglietti congelati, elezioni con scadenza, memorie scritte fuori dalla cellula, misure assolute che rompono gli accordi sbagliati: sono comparse più volte, in forme diverse, sempre per lo stesso motivo.

Un'algebra serve a tre cose:

1. **Scrivere** una forma nuova componendo pezzi già provati, invece di inventarla da capo.
2. **Prevedere** quando una composizione regge e quando no, prima di simularla.
3. **Confrontare** programmi diversi per la stessa forma: due anelli identici a riposo possono avere robustezze opposte.

La domanda di fondo è se gli operatori di composizione conservano la stabilità delle primitive. Se sì, una forma complessa è stabile quanto i suoi pezzi; se no, bisogna sapere esattamente dove si rompe.

## 2. Il modello

Il tessuto è un **grafo dinamico**: i nodi sono le cellule, gli archi i contatti. A differenza dei modelli di calcolo distribuito classici, il grafo cambia mentre il programma gira: le cellule nascono, muoiono, si muovono e quindi cambiano vicine.

**Stato di una cellula.** Ogni cellula *i* porta uno stato

```latex
\sigma_i = (x_i,\; m_i,\; \phi_i,\; F_i,\; M_i)
```

dove x è la posizione, m lo stato della macchina (PROLIFERA, QUIESCENTE, BORDO, …), φ la fase del programma, F i valori dei campi (ρ, d₁, il contatore di quiete, i biglietti) e M le marche (POLO1, DISCO, …).

**Vicinato.** N(i) è l'insieme delle cellule a contatto con *i*. Oltre alle vicine, la cellula può leggere l'ambiente locale: la matrice extracellulare nel volumetto attorno a sé.

**Aggiornamento.** A ogni tick tutte le cellule si aggiornano insieme con la stessa funzione:

```latex
\sigma_i' = \Phi\big(\sigma_i,\ \{\sigma_j : j \in N(i)\},\ E(x_i)\big)
```

dove E(x) è la matrice letta attorno alla cellula. Poi la fisica sposta le cellule e ricalcola i contatti.

**Assioma di località.** Φ non legge nient'altro: né il numero totale di cellule, né la posizione di cellule lontane, né un orologio globale. Ogni informazione globale deve arrivare propagandosi di vicina in vicina, o scritta nella matrice. Le prove hanno trovato due violazioni involontarie di questo assioma (un limite al numero di cellule e un controllo "esiste già un polo?"), ed entrambe nascondevano un problema vero.

## 3. I campi

Un **campo** è una grandezza che ogni cellula conosce di sé e calcola solo dalle vicine. Le prove ne hanno usati tre tipi, con regole di propagazione diverse.

| Tipo | Regola | Esempi | Si stabilizza? |
| --- | --- | --- | --- |
| Distanza | 0 sulla sorgente, altrove il minimo tra le vicine più la distanza da esse | ρ, d₁, d₂ | Sì, su un grafo fisso, in un numero di tick pari al diametro |
| Contatore | 0 dove succede qualcosa, altrove il minimo tra le vicine più uno | quiete | Sì; dice da quanti contatti non c'è attività |
| Biglietto | Il massimo visto tra le vicine, con un'età che cresce a ogni passaggio e scade | elezioni | Sì, se i biglietti sono congelati |

Il campo distanza si scrive:

```latex
f_i = \begin{cases} 0 & i \in S \\ \min_{j \in N(i)} \big(f_j + \lVert x_i - x_j \rVert\big) & \text{altrimenti} \end{cases}
```

Tre fatti su questi campi sono emersi dalle prove e reggono tutto il resto:

- **Se la sorgente sparisce, il campo non si ferma.** I valori salgono senza limite (il "conta all'infinito" dei protocolli di instradamento). Serve un'àncora: un polo che passa la marca a una vicina, una memoria scritta nella matrice.
- **Le differenze si propagano meglio delle medie.** Fare la media di punti stimati che giacciono su un arco li sposta verso l'interno dell'arco: con un errore di direzione δ, ogni passaggio accorcia ρ di circa ρδ²/2, e l'errore si accumula. Propagare ρ per differenze (ρ mio = ρ della vicina meno lo spostamento verso il centro) non ha questo errore sistematico.
- **Un campo che solo si copia non ha àncora.** Se ogni cellula aggiusta ρ solo rispetto alle vicine, l'intero tessuto può slittare di una costante senza accorgersene. Serve almeno una misura assoluta: la curvatura delle fibre, il nucleo, la memoria.

## 4. Le primitive

Ogni primitiva è una tripla di Hoare con in più un evento di fine:

```latex
\{P\}\ p\ \{Q\} \quad \text{con evento di fine } e_p
```

Se il tessuto soddisfa la precondizione P, la regola locale di p lo porta in uno stato che soddisfa la postcondizione Q, e almeno una cellula può rilevare localmente l'evento eₚ che dice "ho finito".

| Primitiva | Precondizione | Postcondizione | Evento di fine |
| --- | --- | --- | --- |
| `campo(c, da S)` | S non vuoto e connesso al resto | ogni cellula conosce la sua distanza da S | il valore non cambia da k tick |
| `cresci(R)` | almeno una cellula in R | R è pieno, nessuna cellula fuori da R | quiete |
| `svuota(R)` | — | nessuna cellula in R | quiete |
| `elegge(P, b)` | b congelato, tessuto connesso | esattamente una cellula ha la marca P | arrivo del campo che parte da P |
| `marca(M, R)` | i campi che definiscono R sono stabili | le cellule in R hanno la marca M | immediato |
| `spingi(M, v)` | M marcato, v definito sui campi | le cellule M si sono spostate lungo v | contatto o quiete |
| `stacca(M, N)` / `attacca(M, N)` | M e N marcati | cambia l'adesione tra M e N | immediato |
| `impegna` | — | nessuna cellula cresce più | immediato |
| `mantieni(R)` | — | R viene ripristinato dopo ogni ferita | mai |

Le primitive non sono tutte dello stesso tipo. `campo`, `marca` ed `elegge` **calcolano**: cambiano campi e marche senza toccare il supporto. `cresci`, `svuota` e `spingi` **costruiscono**: fanno nascere, morire o muovere cellule, quindi cambiano il grafo su cui le altre calcolano. Il secondo gruppo è ciò che distingue questa algebra da quelle del calcolo distribuito su supporto fisso.

## 5. La semantica degli operatori

| Operatore | Significato | Condizione perché regga |
| --- | --- | --- |
| `A ; B` | B parte quando A ha emesso il suo evento di fine | Q\_A implica P\_B, e la fine di A è rilevabile localmente |
| `A ‖ B` | le regole di A e B sono attive insieme | A e B non scrivono le stesse variabili sulle stesse cellule, oppure le loro scritture commutano |
| `A @ M` | A vale solo per le cellule con marca M | M è stabile finché A lavora |
| `A fino a E` | A finché qualche cellula rileva l'evento E | E è rilevabile localmente |
| `sempre A` | A non finisce mai | A non contraddice ciò che fanno gli altri operatori attivi |

**La sequenza come onda.** In un sistema senza orologio globale, "B parte quando A ha finito" non può avvenire ovunque nello stesso istante. La fase φ di ogni cellula avanza in due modi: quando la cellula rileva da sé l'evento di fine, oppure quando tocca una vicina già avanti. La fase è **monotona**, non torna mai indietro, e si diffonde come un'onda alla velocità di un contatto per tick. La regola di sequenza diventa quindi:

```latex
\frac{\{P\}\,A\,\{Q\} \qquad Q \Rightarrow P' \qquad \{P'\}\,B\,\{R\} \qquad e_A \text{ rilevabile}}{\{P\}\ A\,;\,B\ \{R\}}
```

Il punto delicato è la quarta premessa. Se l'evento di fine viene rilevato troppo presto in un settore, l'onda di B parte mentre A lavora ancora altrove. È il problema classico della rilevazione di terminazione nei sistemi distribuiti.

**Il parallelo come non interferenza.** Due regole in parallelo reggono se nessuna invalida ciò che l'altra ha stabilito: è la condizione di non interferenza di Owicki e Gries. L'esempio tipico di interferenza è `sempre mantieni(R) ‖ spingi`: la manutenzione ricostruisce R mentre la spinta lo deforma.

## 6. Le leggi algebriche proposte

Sono congetture: alcune sono quasi ovvie, altre vanno provate nel simulatore o dimostrate.

| Legge | Enunciato | Stato |
| --- | --- | --- |
| Associatività della sequenza | `(A ; B) ; C = A ; (B ; C)` | Attesa, dalla monotonia della fase |
| Commutatività del parallelo | `A ‖ B = B ‖ A`, se A e B non interferiscono | Attesa |
| Idempotenza | `cresci(R) ; cresci(R) = cresci(R)` e lo stesso per `svuota` | Da provare: è un test di "pulizia" della primitiva |
| Manutenzione | `mantieni(R) = sempre (cresci(R) ‖ svuota(¬R))` | È la definizione; regge nelle prove sulle ferite |
| Distributività della restrizione | `(A ; B) @ M = (A @ M) ; (B @ M)`, se M non cambia | Attesa |
| Elemento assorbente | `impegna ; A = impegna` per ogni A che costruisce | Vera per definizione |
| Conflitto | `sempre mantieni(R) ‖ spingi(M, v) = ⊥`, a meno che la spinta tenga le cellule dentro R | Osservata: la manutenzione disfa le deformazioni |
| Chiusura dell'elezione | `elegge(P, b) ; elegge(P, b) = elegge(P, b)` | Vera solo se l'elezione si chiude; senza chiusura le prove davano un secondo vincitore |

La legge di conflitto è la più importante per le trasformazioni. Dice che un tessuto non può essere insieme perfettamente riparabile e deformabile: durante una metamorfosi la manutenzione va sospesa, oppure il suo bersaglio deve spostarsi insieme alla deformazione.

La spinta l'ha confermata in una seconda forma. Un disco che spinge un tessuto pieno non lascia un canale, né se il tessuto è liquido né se è plastico: il materiale spostato gli scorre intorno e riempie la scia. Il volume del canale deve sparire, non spostarsi. Per questo le cellule che il disco trova davanti muoiono (cavitazione), e solo così il foro esce largo: gioco 13–16 contro 5–9 con la sola spinta.

Dopo il foro, ogni cellula stampa la propria posizione nella matrice, e da lì vale di nuovo `mantieni(R) = sempre (cresci(R) ‖ svuota(¬R))`. Cambia solo da dove viene R: non più calcolata da ρ, ma ricordata dalla matrice. Una forma ottenuta per deformazione non ha una formula semplice nelle coordinate della cellula, ma può essere memorizzata. Se una ferita cancella anche lo stampo, quella parte non ricresce: la memoria è l'unica fonte del bersaglio.

## 7. Tre nozioni di robustezza

| Nozione | Domanda | Come si prova |
| --- | --- | --- |
| Auto-stabilizzazione | Da uno stato qualsiasi che rispetta la precondizione, si arriva alla postcondizione e ci si resta? | Crescita da una cellula, molte ripetizioni |
| Resilienza | Dopo una ferita (cellule tolte, matrice cancellata, taglio) si torna alla forma? | Ablazioni, metà tessuto, taglio in due |
| Robustezza al rumore | Con bussole rumorose o cellule che ruotano, la forma regge? | Errore di bussola, rotazione delle cellule |

Le tre nozioni non vanno sempre insieme. Due programmi per l'anello lo mostrano:

- `anello₁ = disco ; svuota(ρ < Rin) ; impegna` è stabile ma non resiliente: dopo `impegna` nessuna cellula ricresce, e una ferita resta aperta.
- `anello₂ = cresci(ρ ≤ R) ‖ svuota(ρ < Rin) @ bloccate` è resiliente: la stessa regola che lo fa crescere lo ripara.

Il compromesso è generale. `impegna` compra stabilità al prezzo della riparazione; `sempre mantieni` compra riparazione al prezzo della deformabilità (la legge di conflitto). Una metamorfosi deve passare per una fase fragile, come fanno gli embrioni veri durante la gastrulazione.

## 8. Regole di buona formazione

Ognuna di queste regole è nata da una prova fallita. Insieme sono i vincoli che un programma deve rispettare perché l'algebra tenga.

1. **Località stretta.** Nessuna regola legge informazioni globali. *Violazione trovata:* un controllo "esiste già un polo?" sull'intero tessuto nascondeva un'elezione che non si chiudeva.
2. **Biglietti congelati.** Un'elezione usa un valore fissato all'ingresso. *Violazione trovata:* con ρ letto a ogni tick, il rumore impediva a qualunque cellula di restare imbattuta.
3. **Elezioni chiuse.** L'arrivo del campo che parte dal vincitore toglie le cellule dalla gara. *Violazione trovata:* le cellule nate dopo per ricambio eleggevano un secondo polo a 130 unità dal primo.
4. **Marche ereditabili.** Una marca su una cellula sopravvive alla cellula, per successione o per memoria nella matrice. *Violazione trovata:* in metà delle prove un polo moriva per ricambio.
5. **Propagazione per differenze.** Le distanze viaggiano come differenze, non come medie di punti. *Violazione trovata:* la media accorciava ρ e il disco esplodeva.
6. **Almeno un'àncora assoluta.** Ogni campo copiato tra vicine ha una fonte che misura direttamente. *Violazione trovata:* senza àncora un intero settore si accordava su un ρ sbagliato.
7. **Fiducia solo dove la lettura è affidabile.** Una cellula legge la matrice solo dove densità e coerenza bastano, e lontano dal centro. *Violazione trovata:* un errore letto con 7 cellule si propagava a tutto il tessuto.
8. **Direzioni mediate tra vicine.** In 3D la direzione del centro viene dalla media dei centri indicati dalle vicine. *Violazione trovata:* gemme, cioè gruppi che si costruivano un secondo centro coerente con se stesso.
9. **Canali che sbagliano in modi diversi.** Le informazioni ridondanti sono utili solo se i loro errori sono indipendenti: dentro la cellula (veloce, deriva), tra vicine (coerente, ma un errore comune non si vede), nella matrice (ferma, ma ricorda gli errori), misura assoluta (rumorosa, ma non copia).

La spinta dei dischi ha aggiunto altre cinque regole, anche queste nate da prove fallite:

10. **Soglia di cedimento.** Un tessuto pieno si comporta come un solido: sotto una certa spinta non cede affatto. *Violazione trovata:* con spinta 0,05 e 0,3 i dischi si appiattivano ma restavano fermi; da circa 0,5 avanzano.
11. **Fine locale, poi onda.** L'evento di fine di una spinta è un contatto che dura nel tempo, e da lì la fase avanza come onda. *Violazioni trovate:* fermare tutto al primo contatto bloccava i dischi a 30 unità, perché un bordo incurvato toccava per primo; non fermarli mai faceva attraversare il tessuto alla coppia di dischi.
12. **Il lato apicale non aderisce.** La coesione forte vale solo tra le cellule già unite al momento dell'`impegna`. *Violazione trovata:* senza questa distinzione le pareti opposte del canale si incollavano e lo chiudevano come una cerniera.
13. **Le cellule intrappolate si tolgono.** Una cellula del corpo che ha i due dischi da parti opposte muore (cucitura). *Violazione trovata:* la regola "tocca entrambi i dischi" uccideva a cascata anche le cellule del bordo, 145 in una prova.
14. **Le marche finiscono.** Una regola utile in una fase può essere dannosa nella successiva. *Violazione trovata:* la successione dei poli trasferiva la marca dal cuscino morente al corpo, e i poli sopravvivevano al foro.

### Invarianza di scala

Un programma morfogenetico ben scritto dovrebbe dare la stessa forma con cellule più grandi o più piccole, come un embrione con più o meno cellule. Per questo ogni parametro va dichiarato in una delle tre unità:

| Unità | Esempi | Come scala |
| --- | --- | --- |
| Lunghezza della forma | R, Rin, raggio del disco, velocità di spinta | fissa |
| Lunghezza della cellula | caselle della matrice, passo di lettura, soglie vicino al centro, rumore, stampo | proporzionale al raggio r della cellula |
| Tempo o contatti | quiete, scadenza dei biglietti, durata delle fasi | fissa, purché resti maggiore del diametro misurato in contatti |

Prova con R = 60 e raggio della cellula 6, 4,5 e 3 (circa 850, 1.900 e 6.100 cellule). Sfera, poli, dischi, foro e riparazione funzionano a tutte e tre le risoluzioni, e l'errore di ρ scende (1,7 → 1,2 → 0,9). Il foro invece non è invariante: in unità assolute si restringe (gioco circa 15 → 12 → 8,5), mentre in raggi di cellula resta circa 3. Qualche regola fissa la larghezza del canale in contatti invece che in lunghezza della forma: è la prossima cosa da trovare.

### Due regole in più, dal blocco della quiete

In 2 prove su 12 la sfera non dichiarava mai la quiete, e tutto il programma restava fermo. La causa era un **ciclo futile**: una cellula appena sotto R si divide, la figlia nasce dentro il bersaglio ma la pressione della superficie la spinge oltre la soglia di morte e muore, e il genitore si ridivide. Un ciclo ogni circa 130 tick, appena sopra la soglia di quiete di 120, con il genitore attivo da 3.700 tick.

15. **Attività è ciò che succede, non ciò che potrebbe succedere.** Il contatore di quiete si azzera su eventi reali (una divisione, una nascita, un rotolamento, l'inizio di una morte), non su uno stato che dichiara una possibilità. *Violazione trovata:* una cellula ferma in PROLIFERA, capace di dividersi ma che non si divideva mai, teneva il tessuto "attivo" per sempre.
16. **Figliata non vitale.** Una cellula la cui figlia muore appena nata smette di dividersi per un po'. *Violazione trovata:* senza, il genitore rifaceva la stessa figlia nello stesso punto all'infinito.

Dopo le due regole: 0 blocchi su 24 prove, a due risoluzioni, e la sfera resta identica (844 cellule, raggio 59–65). Un primo tentativo di correzione — imporre un margine tra "la figlia può nascere" e "la figlia muore" — toglieva il blocco ma rimpiccioliva la sfera del 12% e la rendeva meno regolare: è stato scartato.

### La gastrula: la composizione non basta

Prima prova vera dell'algebra: scrivere una forma nuova cambiando solo il programma, senza aggiungere regole. Il programma è il toro con due argomenti diversi — un polo invece di due, la spinta verso il centro invece che verso l'altro polo, e come evento di fine ρ piccolo invece del contatto:

```
gastrula = sfera ; attendi quiete ;
  P1 = elegge(1, biglietto = ρ) ; campo(d₁, da P1) ;
  marca(D1, d₁ ≤ r_d ∧ superficie) ; impegna ‖ marca(EPI) ;
  ( spingi(D1, verso il centro) ‖ piatto(D1) ‖ svuota(davanti a D1) )
      fino a ρ(D1) ≤ ρ_min ;
  stampa(corpo) ; sempre mantieni(stampo)
```

**Metà riuscita.** La catena gira al primo colpo, 4 prove su 4: il polo viene eletto, il disco marcato, e il disco scende dalla superficie (profondità 56–58) fino al centro (1–15). L'invaginazione avviene. Nessuna regola nuova è servita.

**Metà fallita.** La bocca si richiude: resta una sfera con un disco sepolto dentro, non una coppa. Tre tentativi di chiuderla, tutti dentro il vocabolario esistente, sono falliti in modi istruttivi:

| Tentativo | In termini di algebra | Esito |
| --- | --- | --- |
| L'epitelio segue il disco | `sempre marca(D1, EPI che tocca D1 da dietro)` | Con 2 contatti dilaga (il disco passa da 15 a 81 cellule e deforma la sfera); con 3 non fa nulla |
| Il disco marca il vuoto che lascia | `stampa(canale)` durante la deformazione, non alla fine | Il canale viene marcato e lo svuotamento parte, ma la bocca si richiude lo stesso |
| Il vuoto ha la precedenza sul pieno | `stampa(canale)` prima di `stampa(corpo)` | Nessun cambiamento misurabile |

**Perché.** È una corsa, non una contraddizione logica. La discesa dura circa 100 tick; il corpo si richiude dietro al disco più in fretta di quanto `svuota(canale)` riesca a liberarlo, e quando `stampa(corpo)` fissa la forma, fissa la forma già chiusa.

**Cosa dice al formalismo.** La condizione di non interferenza di Owicki e Gries è logica: due regole reggono se non scrivono le stesse variabili. Qui le regole non si contraddicono affatto, e il programma fallisce lo stesso perché una è più veloce dell'altra. Serve quindi una condizione in più sull'operatore `‖`, di natura quantitativa: **quando due primitive agiscono sulla stessa regione, una che costruisce e una che disfa, il risultato dipende dal rapporto fra le loro velocità.** È la prima legge dell'algebra che non si può enunciare senza numeri, e va aggiunta alla sezione 5 quando sarà formulata bene.

Il toro non ne è toccato: lì il foro non nasce da un canale tenuto aperto, ma dai dischi fusi che muoiono. Verificato che la catena del toro resta identica dopo queste modifiche.

### Il lume: due meccanismi provati, uno solo tiene

Il guscio e l'anello avevano un difetto dichiarato: una **impalcatura**, cioè uno stato di cellula dedicato che teneva aperto il lume e aveva privilegi di crescita suoi (cresce solo verso l'esterno, non riempie i buchi). Un dispositivo arbitrario: non si saprebbe cosa chiedere a una cellula vera per ottenerlo. Si è provato a sostituirlo.

**Primo tentativo: limite di diffusione.** La cavità emerge da sé — una cellula troppo lontana da una superficie bagnata dal mezzo non riceve abbastanza e muore, come il centro necrotico di uno sferoide. Nessun parametro di lume: solo una lunghezza del mezzo. **Non funziona, e il motivo è strutturale: la regola cancella la propria premessa.** Appena il centro muore, il tessuto che lo circondava non è più profondo, ricresce, muore di nuovo. In 6 prove a tre soglie diverse (20, 30, 45) nessuna ha raggiunto la quiete e nessuna ha aperto un buco: solo logoramento perpetuo. Il limite di diffusione spiega come una cavità *si apre*, non come *resta aperta*. (Durante la prova è emerso un errore vero: il campo di profondità era inizializzato a zero su tutte le cellule, e un campo di distanza min-plus inizializzato in basso resta incollato a un salto. Corretto — inizializzazione a MAI e discesa — la profondità misura bene, 0 al bordo e 54 al centro di un disco di raggio 60; ma la conclusione non cambia.)

**Secondo tentativo, quello adottato: la fascia a due bordi.** Il lume è il *fuori-fascia interno*. Una cellula sotto Rin è fuori fascia esattamente come una che supera R, e muore per la stessa ragione. Non è un dispositivo nuovo: è la primitiva della fascia applicata due volte. Serve una sola clausola in più, per l'avvio — un anello deve crescere *attraverso* il proprio lume futuro per arrivare alla fascia — e la clausola è la stessa su entrambi i lati: **finché sono giovane posso stare sotto Rin e mettere una figlia più in fuori; passata la giovinezza il fronte è arrivato, e chi è rimasto indietro è fuori fascia.** Nessuno stato dedicato, nessun privilegio, nessuna cavità dichiarata.

**Terzo risultato, negativo e utile.** Per tenere aperto il lume erano stati aggiunti cinque meccanismi epiteliali: divisione planare (la figlia resta nel foglio, non attraversa la faccia apicale), deposizione basale della lamina (la matrice va sotto, non attorno: il lume resta privo di matrice per sempre), adesione apicale debole (due facce che guardano il lume non si incollano), estrusione per compressione e anoikis, e la distinzione fra una superficie bagnata dal mezzo e un semplice interstizio (misurata: in un tessuto denso l'interno sta a 0.00 di frazione libera, un bordo vero a 0.23–0.34). Spenti uno per uno sull'anello 2D, **nessuno dei cinque cambia il risultato oltre il rumore fra prove.** La fascia fa tutto il lavoro. I cinque restano nel motore solo perché servono alle *deformazioni* (`costringi` ha bisogno della polarità apicale), e lì vanno rimisurati con lo stesso criterio: se non pesano, escono.

Misure finali, un solo motore per 2D e 3D, con ferita (asportato un quarto) e ricrescita:

| forma | n | bordo esterno | bordo interno | errore ρ | quiete | dopo la ferita |
| --- | --- | --- | --- | --- | --- | --- |
| sfera 3D | 841–853 | 59–66 | — | 1.3–1.8 | 3/3 | 3/3 riparate |
| disco 2D | 216–221 | 55–65 | — | 0.8–1.2 | 3/3 | 3/3 riparate |
| guscio 3D | 678–697 | 53–67 | 28–45 | 1.5–2.0 | 3/3 | 3/3 riparate |
| anello 2D | 150–152 | 52–64 | 30–44 | 1.3–2.0 | 3/3 | 3/3 riparate |

Le forme cave restano più rumorose delle piene: il bordo interno è tenuto da una condizione sulla storia della cellula (la giovinezza), non da una misura istantanea, e si vede.

### La gastrulazione senza spingi e senza stampo

Catena, tutta di primitive già in elenco: `cresci(anello) ; elegge(P1) ; campo(d1) ; marca(arco) ; costringi ; sempre mantieni(costringi)`. Il polo è eletto localmente (biglietto congelato = ρ, chiusura). L'arco si allarga mentre la piega scende — l'anello che scendendo si allarga e forza la ripiegatura.

**Tre correzioni trovate per strada, tutte di meccanismo.**

1. **La polarità era quasi sempre spenta.** Il test apicale confrontava la densità di matrice oltre la faccia libera con una soglia assoluta: ma una cellula depone matrice attorno a sé, e la soglia non veniva mai superata. In un anello assestato una cellula su dieci era polarizzata. **Questo invalida l'ablazione del paragrafo precedente**: i cinque meccanismi epiteliali non risultavano inerti, erano semplicemente mai attivi. La polarità è un *confronto fra i due lati*, non una soglia; e in un monostrato nemmeno il baricentro dello spazio libero basta, perché la cellula ha facce libere sui due lati e la media si annulla. La versione che funziona legge il **gradiente della matrice** attorno a sé (8 direzioni nel piano, 14 nello spazio) e punta l'apicale dove la lamina manca: confronto relativo, nessuna unità, e si rinforza da sé perché la deposizione è basale. Risultato: 60 cellule polarizzate su 79 in un monostrato, con l'orientamento giusto (verso fuori sul bordo, verso il lume all'interno). Rifatta l'ablazione con la polarità accesa, sull'anello statico paga solo l'**adesione apicale debole** (senza, 2 prove su 3 non raggiungono la quiete); divisione planare e deposizione basale sono neutre o costano un po' di tempo di quiete.
2. **La crescita dipendente dall'ancoraggio era biologicamente falsa come l'avevo scritta.** Pretendeva matrice *preesistente* dove sarebbe nata la figlia; ma un epitelio che si espande depone lamina nuova man mano. Il costo misurato era un bordo rientrato di una dozzina di unità (48–63 invece di 52–64). Togliendola l'anello torna in misura. L'ancoraggio vincola la **sopravvivenza**, non l'espansione.
3. **La costrizione apicale è sub-cellulare, e fra sfere non si fa un cuneo.** Una cellula che costringe l'apice diventa a cuneo — stretta sopra, larga sotto. Con una sola lunghezza di riposo per contatto, la coppia si accorcia simmetricamente e non produce nessuna coppia di forze. Quindi in un monostrato `costringi` non piega niente: accorcia il foglio, non lo curva. Serve che l'arco costretto sia una **fila distinta** da quella sotto — la striscia bimetallica. Parete di due cellule (Rin 45, R 60) e la piega si fa. Serve anche che l'arco sia **contiguo**: a chiazze la costrizione fa nodi locali e basta. Per la fila esterna l'informazione robusta non è il gradiente di matrice ma la faccia libera che guarda via dal centro virtuale — due cose che la cellula ha già.

**Il risultato, 3 prove su 3.** Il fronte marcato scende dalla superficie (ρ≈57) a ρ 36–42, il canale vuoto lungo l'asse si apre fino a ρ 45–50, il tessuto resta in un pezzo, lo strain delle coppie costrette si stabilizza attorno a −0.24. A KC=0.8 l'arco si autodistrugge; a KC=0.4 metà delle prove non parte.

**E qui si capisce cosa stava coprendo lo `stampo`.** Rilasciata la costrizione, 3 prove su 3 la coppa **si riapre del tutto in 200–800 tick**: la piega è puramente elastica, e la manutenzione di banda si rimangia il solco. Lo stampo non serviva a ricordare la forma: serviva a impedire alla manutenzione di disfarla. La sostituzione onesta non è una memoria ma una **regola mantenuta**: `sempre mantieni(costringi)`. Con l'arco congelato e la costrizione accesa, 3 prove su 3 la coppa sta ferma per 4000 tick senza deriva (fronte 36–42 invariato, pezzi 1).

Una sottigliezza che è costata una lacerazione: `mantieni` tiene accesa la **regola**, non la **marca**. Congelando la marca, le due pareti del solco — entrambe marcate, e ora in contatto — si costringono l'una contro l'altra: le coppie passano da 30 a 199 e il tessuto si spacca in due. Rivalutando il predicato a ogni tick, una cellula che non è più sulla superficie libera smette da sola. È un caso particolare di una regola generale: **una marca è un'ipotesi sul presente, non un ricordo**.

| sostituzione | al posto di | esito |
| --- | --- | --- |
| `costringi` (coppie apicali contigue, fila esterna di una parete a due) | `spingi` (forza imposta verso il centro) | piega 3/3, tessuto integro |
| `sempre mantieni(costringi)` | `stampo` (griglia di voxel che ricorda la forma) | coppa ferma 4000 tick, 3/3 |
| polarità dal gradiente di matrice | soglia assoluta di densità | 60/79 polarizzate invece di 1/10 |
| ancoraggio solo sulla sopravvivenza | ancoraggio anche sulla crescita | bordo torna in misura |

**Domanda aperta che ne nasce.** La coppa regge *finché la regola è accesa*. In un embrione vero la piega a un certo punto si stabilizza da sé, e il candidato è che il tessuto ripiegato **cambi riferimento**: le cellule invaginate rieleggono un proprio centro virtuale, e la loro forma diventa il bersaglio di una banda nuova invece di uno scarto dalla vecchia. Sarebbe `elegge` e `campo` applicate a una sottopopolazione — nessuna primitiva nuova — e permetterebbe di *spegnere* la costrizione. Non provato.

### Spegnere la costrizione: la sottopopolazione differenziata

La domanda era se dare al tessuto invaginato una regola propria permetta di **spegnere** la costrizione invece di tenerla accesa per sempre. **Sì.**

Prima però un esperimento che dice *perché* la coppa si richiude. Rilasciando la costrizione con ogni divisione e ogni morte congelate, 1 prova su 2 la coppa tiene (canale fermo a 43) e l'altra rincula comunque. Quindi le due cause coesistono: **ricrescita** e **rinculo elastico**. Non è solo meccanica, e questo lascia spazio a una regola.

**La versione minima funziona.** Le cellule che hanno costretto diventano una sottopopolazione differenziata: post-mitotiche, e non più giudicate dalla banda della forma madre. La loro regola è "resta come sei" — identità intrinseca, acquisita localmente dal proprio stato, permanente. Con la costrizione **spenta**, 3 prove su 3 la coppa resta ferma 4000 tick (canale 43, 51, 54; dominio di 40, 42, 32 cellule che né si dividono né muoiono). La moria che si vede nei numeri (103 → 78 cellule) è il tessuto *non* differenziato che la piega aveva trascinato sotto Rin: appena la manutenzione riprende muore, ed è corretto — solo il dominio ha il diritto di stare fuori banda.

Questo **sostituisce `impegna`**, che era un interruttore di fase globale acceso durante la metamorfosi e poi spento. Qui non c'è nessuna fase: c'è un tipo cellulare.

**Due alternative provate e scartate.** Il *rimodellamento giunzionale* — un contatto tenuto accorciato abbastanza a lungo adotta la nuova lunghezza come lunghezza di riposo, per contatto e non per cellula — non tiene la coppa: il tessuto si affolla (n da 100 a 113, strain a −0.48) e il solco si richiude, 2 prove su 2. Né serve combinarlo col dominio.

**Il prezzo, misurato.** Asportata metà del dominio, **la coppa non si ripara**: il dominio scende da 36 a 13 cellule e resta lì, il tessuto circostante ricresce e riempie il solco, la coppa è perduta, 2 prove su 2. Il compromesso è ora esplicito e simmetrico:

|  | memoria | riparazione |
| --- | --- | --- |
| `stampo` (griglia di voxel) | globale, nessuna cellula la possiede | sì |
| `sempre mantieni(costringi)` | nessuna, ma il motore resta acceso | sì, finché acceso |
| dominio differenziato, "resta come sei" | nessuna, motore spento | **no** |

**Il campo proprio dovrebbe comprare entrambe**, e l'ho provato: il dominio elegge un polo fra i suoi (biglietto = numero di vicine del dominio), costruisce la propria geodetica dₘ attraverso le sole cellule del dominio, e vive su una fascia \[0, RD\]. Non funziona ancora — la crescita del dominio non riparte dopo la ferita — ma due cose sono state imparate lungo la strada, e valgono più del risultato mancato.

1. **La fascia di un dominio non si può applicare prima che il suo campo sia converguto.** Alla prima iterazione tutte le distanze sono ignote, e "ignota" non è un motivo per morire: nella prima versione il dominio si estingueva in 200 tick. È lo stesso errore che mi era costato il campo di profondità. Vale come regola: *una distanza ancora in discesa non autorizza nessuna decisione irreversibile.*
2. **RD non è un parametro libero.** Dev'essere il raggio a cui l'arco si è fermato — lo stesso numero della fase che l'ha costruito, misurato in un altro riferimento. Un dominio eredita la propria taglia dal processo che l'ha fatto, non la riceve come dato. È un vincolo nuovo sull'algebra: le due primitive in sequenza condividono un parametro, e la composizione non è libera.

### La catena in 3D

Portata sul motore unico, la catena 3D ha fatto emergere tre difetti che in 2D non si vedevano, tutti e tre nelle primitive e non nel port.

1. **Successione del polo.** Con 600 cellule il ricambio uccide P1 più volte (fino a 20 volte in una corsa) e il campo d₁ perde l'ancora: l'arco svanisce a metà della costrizione. Una marca non deve morire con la cellula che la porta — il polo la passa alla vicina viva più prossima.
2. **Chiusura ereditata.** L'elezione non concludeva mai: il ricambio immette continuamente biglietti nuovi che battono il vincitore e nessuno arriva a 45 tick imbattuto. Una cellula nata dove l'elezione è già decisa — cioè con una vicina già chiusa — non si candida.
3. **Astensione di chi non ha vicine.** Una figlia appena nata senza contatti non sa nulla dell'elezione; se si candida col proprio ρ (alto, è in superficie) batte il vincitore. Chi non ha informazione si astiene.

**E un risultato negativo sulla metodologia: la sezione 2D non predice il 3D per la flessione.** Un arco piano è curvo in una direzione, una calotta sferica in due, e resiste alla flessione molto di più. Con i parametri che in 2D davano una coppa profonda, in 3D il fronte non si muove (resta a ρ 52–54, calotta di 12–28 cellule). Il 2D resta utile per trovare i meccanismi — e è stato decisivo per `costringi` — ma non per tararli.

Con costrizione più forte (0.6) e calotta più larga, in 3D l'invaginazione si fa: il fronte scende da ρ≈57 a 38 e poi 31, il canale si apre fino a ρ 43, il tessuto resta in un pezzo, e sotto `mantieni` la coppa sta ferma 4000 tick. **Ma il conto è salato**: le cellule scendono da 643 a 418 e il guscio si restringe (bordo 39–60 invece di 51–67). La calotta, costringendosi, tira dentro tutto il guscio. In 2D non succedeva. Difetto dichiarato, non risolto.

## 12. Le primitive come operazioni sulla simmetria

Una riformulazione, nata da una domanda sul formalismo cristallografico. Si rivela più produttiva di un semplice cambio di linguaggio: dà una legge di composizione, una semplificazione concreta di una primitiva, e due previsioni falsificabili.

**La distinzione da fare prima di tutto.** Non si parla della simmetria della *configurazione* — quella è sempre banale, le cellule stanno in posti precisi e c'è rumore — ma della simmetria della **specifica**: l'insieme delle regole con i loro parametri. Una specifica ha un gruppo G se è invariante sotto G.

Con questa distinzione, ogni primitiva è di uno di due tipi:

- **equivariante**: commuta con l'azione del gruppo. Ruotare i dati iniziali ruota l'esito. Non toglie simmetria, la trasporta.
- **rompente**: seleziona qualcosa che la configurazione non determina.

E le elezioni sono il secondo tipo — ma non sempre, e il *quando* è la legge.

### La legge

> **`elegge(biglietto = f)` rompe la simmetria se e solo se f è costante sulle orbite del gruppo attuale.** Altrimenti è equivariante, e seleziona l'orbita dove f è massima.

Il motivo è immediato: se f è costante sull'orbita, tutti i candidati sono equivalenti e il vincitore è scelto dal rumore. **L'elezione è il dispositivo che converte rumore in scelta** — che è precisamente la rottura spontanea di simmetria della fisica. Se invece f varia sull'orbita, il massimo è determinato e non c'è nessuna scelta da fare.

Applicata alla catena del toro:

| primitiva | gruppo prima → dopo | tipo |
| --- | --- | --- |
| `cresci(fascia R)` | O(3) → O(3) | equivariante: la fascia è invariante per rotazioni |
| `elegge(P1, biglietto = ρ)` | O(3) → C∞v | **rompente**: ρ è costante su tutta la superficie, decide il rumore |
| `campo(d₁)` | C∞v → C∞v | equivariante: d₁ è determinato da P1 |
| `elegge(P2, biglietto = d₁)` | C∞v → C∞v | **equivariante**: d₁ è massimo all'antipodo, un punto solo |
| `marca(disco: d < a)` | C∞v → C∞v | equivariante |
| `costringi(arco)` | C∞v → C∞v | equivariante |

### Tre conseguenze, tutte verificabili

**1. Il secondo polo non è un'elezione.** Il suo biglietto è d₁, che ha un massimo determinato: l'antipodo. Quindi tutto l'apparato dell'elezione — biglietto congelato, età che scade, chiusura, i 45 passi imbattuto — **non serve**. Basta un predicato locale: *sono un massimo locale di d₁ e nessuna vicina ne ha uno più alto*. Resta una scelta residua fra le poche cellule quasi-equivalenti attorno all'antipodo, ma è una scelta locale fra equivalenti, non una rottura di simmetria globale. È una semplificazione concreta di una primitiva, predetta dalla teoria e non trovata per tentativi.

**2. Le due elezioni non possono stare in parallelo.** Previsione falsificabile: eseguendo `elegge(P1) ‖ elegge(P2)` entrambi i biglietti sarebbero ρ, costante sull'orbita, quindi **entrambi i poli deciderebbero dal rumore e non sarebbero antipodali**. L'ordine della catena non è una comodità implementativa: è imposto dalla simmetria. Si prova a costo quasi nullo nel codice attuale, perché le elezioni non dipendono dalla geometria.

**3. Il nostro toro non è esattamente specchiato.** Un toro di rivoluzione ha gruppo D∞h: asse, piani che lo contengono, e **il piano equatoriale**. Quel piano esiste solo se i due poli sono equivalenti. Nella nostra catena non lo sono — P1 è eletto per primo e il suo campo è il biglietto di P2 — quindi l'oggetto realizzato ha C∞v, e la simmetria equatoriale è solo approssimata. Previsione: **le due bocche del foro differiscono in modo sistematico, non casuale.** Nei registri delle prove i due dischi si comportavano effettivamente in modo un po' diverso; va misurato per bene.

### La legge di conteggio

Se tutte le primitive sono equivarianti sotto G e la specifica iniziale è G-invariante, allora per induzione ogni specifica successiva è G-invariante: nessuna catena di sole primitive equivarianti può specificare un bersaglio meno simmetrico del suo punto di partenza. Quindi:

> **Il numero di elezioni rompenti necessarie in una catena è almeno la profondità della discesa dal gruppo iniziale a quello del bersaglio.**

E i conti tornano, con una sorpresa: disco e sfera **zero** elezioni; gastrula (C∞v) **una**; toro (D∞h realizzato come C∞v) **una**, non due. Un corpo con due assi indipendenti — una simmetria bilaterale — ne richiederebbe **due**, e la seconda dovrebbe usare un biglietto costante sul cerchio lasciato libero dalla prima: eleggere un punto sull'equatore.

Questo è un candidato serio per la legge che mancava — quella che altrove in questo quaderno è annotata come "la legge quantitativa del `‖`". Non è ancora quella: dice *quando* due operazioni non possono essere parallele, non *quanto* devono distare quando possono. Ma è il primo enunciato dell'algebra che vieta qualcosa per ragioni di principio invece che per esperienza.

**Una lettura in più.** Anche il dominio differenziato si dice in questo linguaggio: è una simmetria rotta **resa permanente localmente**, senza bisogno che un campo la ricordi. La differenziazione è il modo in cui una rottura di simmetria diventa stabile.

### Le prove

Misurato sul motore a sfere (2D, R=60, r=4, 8 ripetizioni per modo): l'angolo fra i due poli, che la legge predice \~180° in sequenza e senza motivo di essere 180° in parallelo.

| modo | biglietto di P2 | angolo fra i poli |
| --- | --- | --- |
| **sequenza** | d₁ (varia sull'orbita) | **170° ± 6°** |
| **parallelo, biglietto identico** | ρ, lo stesso per entrambe | **0° in 8 prove su 8**: vince la stessa cellula |
| **parallelo, errore indipendente** | ρ + errore di misura proprio | **61–80°, dev ≈ 50–58°**: sparsi |

Le prime due conseguenze sono confermate, e la seconda più di quanto fosse stato predetto: `elegge(P1,ρ) ‖ elegge(P2,ρ)` non dà due poli non antipodali, **non dà due poli affatto** — con biglietto identico la stessa cellula vince entrambe le elezioni. Il parallelo è degenere, non solo sbagliato.

E la versione con errore indipendente ha un **passaggio al limite misurato**. La variazione residua di ρ sulle cellule di superficie è di 3.2 unità di deviazione (estremi 50–61). Con errore 0.3, molto minore, 5 prove su 8 eleggono ancora la stessa cellula (media 38°): entrambe le elezioni cadono sulla stessa gobba. Con errore 8 e 25, molto maggiori, i due vincitori sono indipendenti e uniformi (medie 80° e 72°, contro i 90° ± 52° attesi per due punti casuali). Quindi: **il parallelo di due elezioni con lo stesso biglietto è deciso da chi è più grande fra la struttura residua del biglietto e l'errore di misura** — e in nessuno dei due regimi si ottengono poli antipodali.

**La terza conseguenza è smentita, e la smentita vale più della previsione.** Avevo dedotto che, essendo la selezione di P2 equivariante, tutto l'apparato dell'elezione fosse inutile e bastasse il predicato locale *sono un massimo locale di d₁*. Provato: conclude più in fretta (122 passi contro 166) ma **si dichiarano da 3 a 8 candidati per prova**, e l'angolo peggiora a 135° ± 34°. In un campo discreto e rumoroso i massimi locali spuri sono molti, e "massimo locale" non è un predicato sufficiente.

Il motivo è che **un'elezione fa due lavori distinti, e la legge di simmetria parla solo del primo**:

1. **rompere la simmetria** quando il biglietto è costante sulle orbite — cioè convertire rumore in scelta;
2. **trovare il massimo globale** di un biglietto che varia ma è rumoroso e ha massimi locali spuri.

Per P2 il lavoro 1 non serve, ma il lavoro 2 sì, e l'apparato (massimo che viaggia, età che scade, chiusura, attesa imbattuto) è lì per quello. La legge resta valida come enunciato sulla simmetria; la semplificazione che ne avevo dedotto no, perché confondeva le due funzioni. Va annotato come distinzione permanente fra le primitive: **rompere ≠ massimizzare.**

Resta in piedi, e con un primo numero, la previsione sull'asimmetria: la sequenza dà 170° ± 6° su 7 prove, cioè un difetto di \~10° rispetto ai 180°, significativo rispetto all'errore della media. L'oggetto realizzato non ha il piano equatoriale: è C∞v, non D∞h, e il difetto è misurabile. Da verificare se si ritrova nelle due bocche del foro del toro.

## 9. Le forme come programmi

```
disco    = campo(ρ, da centro) ; cresci(ρ ≤ R) ; sempre mantieni(ρ ≤ R)

anello₁  = disco fino a quiete ; svuota(ρ < Rin) ; impegna
anello₂  = campo(ρ) ; (cresci(ρ ≤ R) ‖ svuota(ρ < Rin) @ bloccate)

ellisse  = campo(asse) ; campo(ρ con la metrica dell'asse) ; cresci(ρ ≤ R) ; sempre mantieni(ρ ≤ R)

sfera    = disco in tre dimensioni
guscio   = anello₂ in tre dimensioni

planaria = disco ; [taglio] ; sempre mantieni(ρ ≤ R) @ ogni frammento
```

La trasformazione sfera → toro, nella parte già provata e in quella da costruire:

```
sfera ;
attendi quiete ;                                                 -- provata
P1 = elegge(1, biglietto = ρ) ; campo(d₁, da P1) ;              -- provata: un solo vincitore
P2 = elegge(1, biglietto = d₁) ; campo(d₂, da P2) ;             -- provata: angolo 155–174°
marca(D1, d₁ ≤ r_d ∧ superficie) ‖ marca(D2, d₂ ≤ r_d ∧ superficie) ;
impegna ‖ marca(EPI, superficie, con giunzioni) ;               -- provata
( spingi(D1, −∇d₂) ‖ spingi(D2, −∇d₁) ‖ piatto(D)
  ‖ svuota(davanti a D) ) fino a contatto stabile(D1, D2) ;    -- provata: il disco scava il canale (cavitazione)
attacca(D1, D2) ‖ cucitura ;                                     -- provata
stacca(D1 ∪ D2, corpo) ;
fora = svuota(D1 ∪ D2) ;                                        -- provata
stampa(corpo) ; sempre mantieni(stampo)                          -- provata: le ferite guariscono
```

La linea tra i due poli è l'asse del toro: la rottura di simmetria nasce dalle due elezioni, senza bisogno di una polarità ereditata.

## 10. Enunciati da verificare

| Livello | Strumento | Enunciato candidato |
| --- | --- | --- |
| Campi su grafo fisso | Lean | Il campo distanza converge da qualunque stato iniziale in un numero di tick non superiore al diametro del grafo |
| Campi su grafo fisso | Lean | Se la sorgente sparisce, i valori crescono senza limite: nessun punto fisso finito |
| Elezione | TLA+ | Con biglietti congelati e distinti, scadenza maggiore del diametro e chiusura, alla fine esiste esattamente un vincitore |
| Elezione | TLA+ | Senza chiusura, l'ingresso di una cellula con biglietto maggiore produce un secondo vincitore (controesempio atteso) |
| Quiete | TLA+ | La quiete viene dichiarata solo se nessuna cellula entro la soglia di contatti è attiva da almeno soglia tick |
| Sequenza | Lean | La regola di sequenza della sezione 5 è corretta per fasi monotone e eventi rilevabili |
| Parallelo | Lean | Due primitive che non scrivono le stesse variabili commutano |
| Fisico | Simulazione | Ogni primitiva e ogni catena regge su molte ripetizioni, con e senza rumore |

La parte fisica (posizioni continue, nascite e morti, forze) resta affidata alla simulazione. Le dimostrazioni conviene farle quando il linguaggio smette di cambiare: ogni definizione rivista costringe a rifare le prove formali.

## 11. Lavori vicini e domande aperte

Verificati il 26 ottobre 2026. Il campo esiste e si chiama **morphogenetic engineering** (Doursat, Sayama, Michel, dal 2009): auto-organizzazione più architettura, cioè ottenere una forma voluta da regole locali.

- **Field calculus e aggregate programming** (Beal, Viroli): il parente più stretto. Poche primitive su campi distribuiti, con teoremi di composizione per i blocchi auto-stabilizzanti. Differenza decisiva: il supporto è fisso, i nodi non nascono, non muoiono e non si spostano. Tutte le loro primitive sono di quelle che *calcolano*; nessuna *costruisce*.
- **Amorphous computing** (MIT): Growing Point Language di Coore, Origami Shape Language di Nagpal, poi Proto. La rassegna di riferimento, *Organizing the Aggregate* (Beal et al.), dice due cose che ci riguardano: i linguaggi di pattern «generalmente non si compongono in modo pulito», e la rassegna non presenta né teoremi di composizione né leggi algebriche.
- **Morphogenetic Primitives** (Bai, Eyiyurekli, Breen, Drexel, 2008–11): nome vicino, cosa diversa. Lì una primitiva è una *particella* (una cellula virtuale), non un'operazione, e "shape composition" significa aggregazione di particelle, non composizione di comportamenti: nessun operatore, nessuna legge. Le regole sono evolute con programmazione genetica da un bersaglio disegnato a mano, invece che scritte. Divisione, morte, adesione e differenziamento sono esplicitamente disattivati, con numero di agenti fisso. Solo 2D; l'unico bersaglio con un buco (un anello) non converge, e l'affidabilità crolla con la complessità (ellisse 100%, ingranaggio 34%, boomerang 4%). I loro lavori futuri dichiarati — divisione, morte, 3D, auto-riparazione — sono il nostro punto di partenza.
- **Programmazione differenziabile per cluster di cellule** (Nature Computational Science, 2025): impara le regole locali da una forma voluta, con morfogeni, adesione e stress meccanico. Approccio inverso al nostro. Non tratta cambi di topologia né composizione di comportamenti.
- **Logica dei programmi concorrenti:** le triple di Hoare e la non interferenza di Owicki e Gries per il parallelo; la rilevazione di terminazione per la sequenza.
- **Biologia dello sviluppo:** informazione posizionale (Wolpert); *target morphology* e memoria anatomica nella planaria (Levin), parente diretto del nostro stampo nella matrice; reciprocità tra cellula e matrice (Bissell).

Lo spazio che resta libero, per quanto si vede da qui: primitive che **cambiano il supporto** (nascita, morte, movimento), composte da **operatori con leggi**, capaci di **cambiare la topologia** (il genere della superficie), con il bersaglio **ricordato fuori dalla cellula** invece che calcolato.

Domande aperte:

- [ ] Le leggi della sezione 6 reggono nel simulatore, una per una?
- [ ] Esiste un insieme minimo di primitive da cui si ottengono tutte le altre?
- [ ] Come si esprime una forma senza simmetria, che chiede due assi o una chiralità?
- [ ] `spingi` e `stacca`/`attacca` si lasciano descrivere come triple di Hoare, o la meccanica richiede una semantica diversa?
- [ ] Si può dare al bersaglio una forma generale (una funzione nelle coordinate locali) così che una metamorfosi sia solo un bersaglio che si sposta nel tempo?
