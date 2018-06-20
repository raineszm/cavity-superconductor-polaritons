import csv
import itertools as it
from pathlib import Path

import numpy as np
import tqdm

from .. import _bardasis_schrieffer as bsm
from ..materials import NIOBIUM
from ..polariton import build_polariton


def vc(params):
    sys0 = bsm.System(vs=0., **params.material)
    state0 = bsm.State.solve(sys0, params.Trel * params.material["Tc"])
    return (
        state0.delta / sys0.kf
    )  # The critical superfluid velocity via Landau's argument


def qc(p: bsm.Polariton):
    return 2 * p.state.delta / bsm.C


def data(fname, qs, thetas, params, hamiltonian=None):
    if hamiltonian is None:
        hamiltonian = params.cls == bsm.ModePolariton

    p = build_polariton(**params.as_args(vc(params)))

    qs *= p.state.delta / bsm.C

    with open(fname, "w") as f:
        writer = csv.DictWriter(f, ["q", "theta", "omega", "i"])
        writer.writeheader()
        for (theta, q) in tqdm.tqdm(
            it.product(thetas, qs), total=len(thetas) * len(qs)
        ):
            if hamiltonian:
                H = p.hamiltonian(q, theta)
                modes = np.linalg.eigvalsh(H)
            else:
                modes = p.find_modes(q, theta)
            for (i, m) in enumerate(sorted(modes)):
                writer.writerow({"q": q, "theta": theta, "omega": m, "i": i})
