import csv
import itertools as it
from pathlib import Path

import numpy as np
import tqdm

from . import _bardasis_schrieffer as bsm
from .materials import NIOBIUM
from .polariton import niobium

TREL = 0.4  # T/T_c
VREL = 0.6  # vs/v_c
DIPOLEX = 100  # Phenominological enhancement factor
ROOTREL = 0.95


def vc():
    sys0 = bsm.System(vs=0., **NIOBIUM)
    state0 = bsm.State.solve(sys0, TREL * NIOBIUM["Tc"])
    return (
        state0.delta / sys0.kf
    )  # The critical superfluid velocity via Landau's argument


def qc(p: bsm.Polariton):
    return 2 * p.state.delta / bsm.C


def data(fname, qs, thetas):
    p = niobium(
        VREL * vc(), Trel=TREL, dipoleX=DIPOLEX, root_rel=ROOTREL, cls=bsm.ModePolariton
    )

    with open(fname, "w") as f:
        writer = csv.DictWriter(f, ["q", "theta", "omega", "i"])
        writer.writeheader()
        for (theta, q) in tqdm.tqdm_notebook(
            it.product(thetas, qs), total=len(thetas) * len(qs)
        ):
            qx = q * np.cos(theta)
            qy = q * np.sin(theta)

            H = p.hamiltonian(qx, qy)
            modes = np.linalg.eigvalsh(H)
            for (i, m) in enumerate(sorted(modes)):
                writer.writerow({"q": q, "theta": theta, "omega": m, "i": i})
