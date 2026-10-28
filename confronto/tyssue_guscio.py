"""Confronto esterno 3D: il guscio a due superfici minimizzato con tyssue (prova F7).

La geometria di partenza e le tensioni, già tradotte nelle convenzioni di tyssue, vengono dal
nostro motore (confronto/esporta_guscio.jl, file guscio_<scenario>_*.csv). Qui si costruisce un
Monolayer di tyssue con quelle facce, lo si porta verso il minimo con il suo risolutore
quasi-statico (L-BFGS-B) tenendo fermi i vertici del bordo (is_active = 0), e si scrivono le
posizioni finali, il gradiente di tyssue in quelle posizioni e l'energia in
guscio_<scenario>_tyssue_vertici.csv. La batteria Julia (test/fisica_confronto_3d.jl) ricalcola
energia e forze con il nostro motore sulle stesse posizioni.

Effettori di tyssue usati, con le loro formule (verificate sul codice, versione 1.1.0):
- SurfaceTension: surface_tension · area, per copia di faccia;
- LineTension: line_tension · lunghezza / 2, per mezzo-spigolo;
- CellVolumeElasticity: ½ · vol_elasticity · (vol − prefered_vol)².
Area e volume: a ventaglio dal baricentro dei vertici della faccia, come il nostro motore.

Si lancia con:  python -I confronto/tyssue_guscio.py   (vedi confronto/LEGGIMI.md)
"""

import csv
import math
import os
import sys
import types

sys.modules.setdefault("tyssue._mesh_generation", types.SimpleNamespace(make_spherical=None))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from tyssue import Monolayer, config  # noqa: E402
from tyssue.dynamics import effectors, model_factory  # noqa: E402
from tyssue.geometry.bulk_geometry import BulkGeometry  # noqa: E402
from tyssue.solvers.quasistatic import QSSolver  # noqa: E402

CARTELLA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "risultati")
SCENARI = ("anello_centrale", "placca")
SEGMENTO = {"apicale": "apical", "basale": "basal", "laterale": "lateral"}

modello = model_factory(
    [effectors.SurfaceTension, effectors.LineTension, effectors.CellVolumeElasticity]
)


def leggi(nome, parte):
    with open(os.path.join(CARTELLA, f"guscio_{nome}_{parte}.csv")) as fh:
        return list(csv.DictReader(fh))


def costruisci(nome):
    vert = leggi(nome, "vertici")
    facce = leggi(nome, "facce")
    cell = leggi(nome, "cellule")
    # tenuti da parte: la costruzione del Monolayer riempie le tabelle che riceve con i suoi
    # valori predefiniti (is_active = 1 compreso), quindi li si rimette dopo da qui
    attivi = [1 - int(r["fisso"]) for r in vert]
    vert_df = pd.DataFrame(
        {
            "x": [float(r["x"]) for r in vert],
            "y": [float(r["y"]) for r in vert],
            "z": [float(r["z"]) for r in vert],
            "is_active": [1 - int(r["fisso"]) for r in vert],
            "segment": ["apical" if int(r["vertice"]) % 2 == 1 else "basal" for r in vert],
        },
        index=pd.Index(range(len(vert)), name="vert"),
    )
    righe_spigoli, righe_facce = [], []
    for fi, r in enumerate(facce):
        v = [int(k) - 1 for k in r["vertici"].split()]
        lam = [float(k) for k in r["tensioni_linea"].split()]
        c = int(r["cellula"]) - 1
        seg = SEGMENTO[r["etichetta"]]
        righe_facce.append({"segment": seg, "surface_tension": float(r["tensione"]),
                            "is_alive": 1, "is_active": 1})
        for k in range(len(v)):
            righe_spigoli.append({"srce": v[k], "trgt": v[(k + 1) % len(v)], "face": fi,
                                  "cell": c, "segment": seg, "line_tension": lam[k],
                                  "is_active": 1})
    edge_df = pd.DataFrame(righe_spigoli)
    edge_df.index.name = "edge"
    face_df = pd.DataFrame(righe_facce)
    face_df.index.name = "face"
    cell_df = pd.DataFrame(
        {"prefered_vol": [float(r["V0"]) for r in cell],
         "vol_elasticity": [float(r["K"]) for r in cell], "is_alive": 1},
        index=pd.Index(range(len(cell)), name="cell"),
    )
    specs = config.geometry.bulk_spec()
    eptm = Monolayer(nome, {"vert": vert_df, "edge": edge_df, "face": face_df, "cell": cell_df},
                     specs, coords=["x", "y", "z"])
    # i valori di tensione e volume vanno rimessi dopo la costruzione, che riempie i default
    eptm.face_df["surface_tension"] = face_df["surface_tension"].values
    eptm.edge_df["line_tension"] = edge_df["line_tension"].values
    eptm.cell_df["prefered_vol"] = cell_df["prefered_vol"].values
    eptm.cell_df["vol_elasticity"] = cell_df["vol_elasticity"].values
    eptm.vert_df["is_active"] = attivi
    BulkGeometry.update_all(eptm)
    assert int((eptm.vert_df["is_active"] == 0).sum()) == attivi.count(0), "bordo non bloccato"
    return eptm


def forza_massima(eptm):
    g = modello.compute_gradient(eptm)
    attivi = eptm.vert_df["is_active"].astype(bool)
    return float(np.sqrt((g[attivi.values] ** 2).sum(axis=1)).max())


def minimizza(eptm, forza_obiettivo=1e-7, giri_massimi=300):
    risolutore = QSSolver()
    for giro in range(1, giri_massimi + 1):
        risolutore.find_energy_min(
            eptm, BulkGeometry, modello,
            options={"gtol": 1e-14, "ftol": 0.0, "maxiter": 200000, "maxfun": 400000},
        )
        BulkGeometry.update_all(eptm)
        f = forza_massima(eptm)
        if f < forza_obiettivo:
            break
    return giro, f


if __name__ == "__main__":
    from importlib.metadata import version
    print("tyssue", version("tyssue"))
    for nome in SCENARI:
        eptm = costruisci(nome)
        print(nome, "energia iniziale:", float(modello.compute_energy(eptm)))
        giri, f = minimizza(eptm)
        E = float(modello.compute_energy(eptm))
        print(nome, "giri:", giri, "| forza massima residua:", f, "| energia:", E)
        scrivi_stato(nome, eptm)


def scrivi_stato(nome, eptm):
    """Posizioni finali e gradiente di tyssue in quelle posizioni (zero sui vertici bloccati):
    la batteria ricalcola energia e forze con il nostro motore sulle stesse posizioni."""
    g = modello.compute_gradient(eptm)
    g = g.reindex(eptm.vert_df.index).fillna(0.0)
    g[eptm.vert_df["is_active"].values == 0] = 0.0
    tab = eptm.vert_df[["x", "y", "z"]].copy()
    tab[["gx", "gy", "gz"]] = g[["gx", "gy", "gz"]].values
    tab["energia"] = float(modello.compute_energy(eptm))
    tab.to_csv(os.path.join(CARTELLA, f"guscio_{nome}_tyssue_vertici.csv"),
               header=True, index=False, float_format="%.17g")
