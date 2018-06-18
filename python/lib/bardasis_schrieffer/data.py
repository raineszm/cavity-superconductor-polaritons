from . import _bardasis_schrieffer as bsm
import csv
import itertools as it
import numpy as np
import tqdm
from pathlib import Path
from polariton import niobium
from materials import NIOBIUM


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
    writer = csv.DictWriter(fname, ["q", "theta", "omega"])
    writer.writeheader()
    for (theta, q) in tqdm.tqdm(it.product(thetas, qs)):
        qx = q * np.cos(theta)
        qy = q * np.sin(theta)

        modes = p.find_modes(qx, qy)
        for m in modes:
            writer.writerow({"q": q, "theta": theta, "omega": m})

