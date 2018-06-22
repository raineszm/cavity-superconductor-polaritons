import csv
import itertools as it
from pathlib import Path
from multiprocessing import Pool

import numpy as np
import tqdm

from .. import _bardasis_schrieffer as bsm
from ..materials import NIOBIUM
from ..polariton import build_polariton


class Runner:
    def __init__(self, p, hamiltonian):
        self.p = p
        self.hamiltonian = hamiltonian

    def __call__(self, args):
        (theta, q) = args
        print(q, theta)
        if self.hamiltonian:
            modes = self.p.bands(q, theta)
        else:
            modes = self.p.find_modes(q, theta)
        return [
            {"q": q, "theta": theta, "omega": m, "i": i}
            for (i, m) in enumerate(sorted(modes))
        ]


def data(fname, qs, thetas, params, hamiltonian=None):
    if hamiltonian is None:
        hamiltonian = params.cls == bsm.ModePolariton

    p = build_polariton(**params.as_args())

    qs *= p.state.delta / bsm.C

    pool = Pool()
    runner = Runner(p, hamiltonian)
    with open(fname, "w") as f:
        writer = csv.DictWriter(f, ["q", "theta", "omega", "i"])
        writer.writeheader()
        for rows in tqdm.tqdm(
            pool.imap_unordered(runner, it.product(thetas, qs)),
            total=len(thetas) * len(qs),
        ):
            writer.writerows(rows)
