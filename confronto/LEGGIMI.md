# Confronto esterno

Lo stesso esperimento in un motore maturo, per validare la fisica (specifica, sezione 9).

## tyssue — modello a vertici laterale (prova F5)

`tyssue_laterale.py` costruisce gli stessi tessuti della prova F5, li porta al minimo
dell'energia con tyssue e scrive le misure per cellula in `risultati/*.csv`. La batteria Julia
(`test/fisica_confronto.jl`) legge quei file: per lanciare le prove non serve Python.

Risultato al 28 ottobre 2026: i due motori coincidono a meno di 1e-6 su angoli, altezze e aree,
e l'energia totale fino a 1e-13 relativo, su tre casi: striscia piatta con una e con quattro
cellule costrette, anello con quattro cellule costrette.

### Rigenerare i risultati

tyssue 1.1.0 non funziona con NumPy 2 (usa il prodotto vettoriale in 2D, che NumPy 2 ha tolto),
e NumPy 1.26 non c'è per Python 3.13. Serve Python 3.12:

```
python3.12 -m venv venv-tyssue
venv-tyssue/bin/pip install "numpy<2" "pandas<2.3" scipy scikit-image quantities tyssue==1.1.0
venv-tyssue/bin/python -I confronto/tyssue_laterale.py
```

Versioni usate: numpy 1.26.4, pandas 2.2.3, scipy 1.17.1, tyssue 1.1.0. Ci vogliono circa sei
minuti: il minimizzatore di tyssue (L-BFGS-B) si ferma da solo con forze residue di 0.01 nN, e
lo script lo fa ripartire finché scendono sotto 1e-7 nN.

Il pacchetto di PyPI non porta l'estensione compilata con CGAL per le forme sferiche: lo script
la sostituisce con un segnaposto, perché qui non serve.

Licenza di tyssue: GPL-3.0, verificata sul pacchetto installato. È usato come programma
separato; il suo codice non entra nel motore.
