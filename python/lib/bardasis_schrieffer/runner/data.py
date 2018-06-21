import csv
import itertools as it
from pathlib import Path

import numpy as np
import tqdm

from .. import _bardasis_schrieffer as bsm
from ..materials import NIOBIUM
from ..polariton import build_polariton


def data(fname, qs, thetas, params, hamiltonian=None):
    if hamiltonian is None:
        hamiltonian = params.cls == bsm.ModePolariton

    p = build_polariton(**params.as_args())

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
