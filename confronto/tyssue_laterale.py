"""Confronto esterno: lo stesso modello a vertici laterale minimizzato con tyssue.

tyssue (DamCB, GPL-3.0, https://tyssue.readthedocs.io) è usato come programma separato: il suo
codice non entra nel nostro motore. Questo script costruisce gli stessi tessuti delle prove,
li porta al minimo dell'energia con il risolutore quasi-statico di tyssue (L-BFGS-B di scipy)
e scrive le misure per cellula in confronto/risultati/*.csv. La prova F5 della batteria Julia
confronta il nostro motore con questi file.

Convenzioni di tyssue da tradurre (verificate sul suo codice, versione 1.1.0):
- LineTension vale line_tension · lunghezza / 2 per mezzo-spigolo. Uno spigolo laterale
  condiviso è fatto di due mezzi-spigoli, quindi line_tension = Γl; uno spigolo apicale,
  basale o laterale di bordo è un solo mezzo-spigolo, quindi line_tension = 2Γ.
- FaceAreaElasticity vale ½ · area_elasticity · (area − prefered_area)², come la nostra K.

Si lancia con:  python -I confronto/tyssue_laterale.py   (con tyssue installato)
"""

import csv
import math
import os
import sys
import types

# La generazione sferica di tyssue richiede un'estensione compilata con CGAL che il pacchetto
# da PyPI non porta. Qui non serve: la si sostituisce con un segnaposto prima dell'import.
sys.modules.setdefault("tyssue._mesh_generation", types.SimpleNamespace(make_spherical=None))

import tyssue  # noqa: E402
from tyssue.dynamics import effectors, model_factory  # noqa: E402
from tyssue.generation.shapes import generate_lateral_tissue, generate_ring  # noqa: E402
from tyssue.geometry.planar_geometry import PlanarGeometry  # noqa: E402
from tyssue.solvers.quasistatic import QSSolver  # noqa: E402

# Parametri condivisi con test/fisica_confronto.jl
GA_COSTRETTA, GA, GB, GL = 1.3, 1.0, 1.0, 1.0
A0, K = 100.0, 0.01
N = 20

CARTELLA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "risultati")
modello = model_factory([effectors.LineTension, effectors.FaceAreaElasticity])


def semiangolo(u, v):
    croce = u[0] * v[1] - u[1] * v[0]
    return math.atan2(croce, u[0] * v[0] + u[1] * v[1]) / 2


def imposta(eptm, costrette, bordo_laterale):
    """Tensioni e area. `bordo_laterale`: indici dei mezzi-spigoli laterali senza gemello."""
    e = eptm.edge_df
    e["is_active"] = 1
    for i in e.index:
        f = int(e.at[i, "face"])
        seg = e.at[i, "segment"]
        if seg == "apical":
            e.at[i, "line_tension"] = 2 * (GA_COSTRETTA if f in costrette else GA)
        elif seg == "basal":
            e.at[i, "line_tension"] = 2 * GB
        else:
            e.at[i, "line_tension"] = 2 * GL if i in bordo_laterale else GL
    eptm.face_df["area_elasticity"] = K
    eptm.face_df["prefered_area"] = A0
    eptm.face_df["is_alive"] = 1


def minimizza(eptm, forza_obiettivo=1e-9, giri_massimi=200):
    """L-BFGS-B di tyssue, fatto ripartire finché la forza massima su un vertice scende sotto
    `forza_obiettivo` (nN): da solo si ferma quando l'energia smette di scendere in modo
    apprezzabile, con forze residue dell'ordine di 0.01 nN."""
    risolutore = QSSolver()
    for giro in range(1, giri_massimi + 1):
        PlanarGeometry.update_all(eptm)
        res = risolutore.find_energy_min(
            eptm, PlanarGeometry, modello,
            options={"gtol": 1e-14, "ftol": 0.0, "maxiter": 200000, "maxfun": 400000},
        )
        PlanarGeometry.update_all(eptm)
        g = modello.compute_gradient(eptm)
        res.forza_massima = float((g ** 2).sum(axis=1).max() ** 0.5)
        res.giri = giro
        if res.forza_massima < forza_obiettivo:
            break
    return res


def scrivi(nome, eptm, celle, energia):
    with open(os.path.join(CARTELLA, nome + ".csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cellula", "semiangolo", "altezza", "area", "energia_totale"])
        for c, (u, v, area) in enumerate(celle, start=1):
            altezza = (math.hypot(*u) + math.hypot(*v)) / 2
            w.writerow([c, repr(semiangolo(u, v)), repr(altezza), repr(area), repr(energia)])


def pos(eptm, k):
    return eptm.vert_df.at[k, "x"], eptm.vert_df.at[k, "y"]


def diff(a, b):
    return a[0] - b[0], a[1] - b[1]


def striscia(nome, costrette):
    # basale y = 0 (vertici 0..N), apicale y = 14 (vertici N+1..2N+1), larghezza 7
    eptm = generate_lateral_tissue(N, 7.0 * N, 14.0)
    lat = eptm.edge_df[eptm.edge_df["segment"] == "lateral"]
    bordo = set()
    for i in lat.index:
        s, t = int(lat.at[i, "srce"]), int(lat.at[i, "trgt"])
        if {s, t} == {0, N + 1} or {s, t} == {N, 2 * N + 1}:
            bordo.add(i)
    imposta(eptm, costrette, bordo)
    res = minimizza(eptm)
    celle = []
    for f in range(N):
        b0, b1 = pos(eptm, f), pos(eptm, f + 1)
        a0, a1 = pos(eptm, f + N + 1), pos(eptm, f + N + 2)
        celle.append((diff(a0, b0), diff(a1, b1), eptm.face_df.at[f, "area"]))
    scrivi(nome, eptm, celle, float(modello.compute_energy(eptm)))
    return res


def anello(nome, costrette):
    # apicale all'interno: vertici 0..N-1 apicali, N..2N-1 basali; tyssue gira in senso orario
    eptm = generate_ring(N, 25.0, 39.0, apical="in")
    imposta(eptm, costrette, set())
    res = minimizza(eptm)
    celle = []
    for f in range(N):
        g = (f + 1) % N
        # nel nostro ordine antiorario il lato "sinistro" della faccia f è quello del vertice g
        u = diff(pos(eptm, g), pos(eptm, g + N))
        v = diff(pos(eptm, f), pos(eptm, f + N))
        celle.append((u, v, eptm.face_df.at[f, "area"]))
    scrivi(nome, eptm, celle, float(modello.compute_energy(eptm)))
    return res


if __name__ == "__main__":
    os.makedirs(CARTELLA, exist_ok=True)
    from importlib.metadata import version
    print("tyssue", version("tyssue"))
    # facce numerate da 0; le costrette sono al centro
    for nome, f, cost in [
        ("striscia_k1", striscia, {9}),
        ("striscia_k4", striscia, {8, 9, 10, 11}),
        ("anello_k4", anello, {0, 1, 2, 3}),
    ]:
        r = f(nome, cost)
        print(nome, "giri:", r.giri,
              "| forza massima residua:", r.forza_massima)
