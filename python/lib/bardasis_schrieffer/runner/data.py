import csv
import itertools as it
from multiprocessing import Pool, Value
from pathlib import Path

import click
import numpy as np
import tqdm

from .. import _bardasis_schrieffer as bsm
from ..materials import NIOBIUM
from ..polariton import build_polariton


class Runner:
    FAILURES = Value("i", 0)

    def __init__(self, p, hamiltonian, params, Nfail):
        self.p = p
        self.hamiltonian = hamiltonian
        self.params = params
        self.Nfail = Nfail

    def __call__(self, args):
        (theta, q) = args
        try:
            if self.hamiltonian:
                modes = self.p.bands(q, theta)
            else:
                modes = self.p.find_modes(
                    q,
                    theta,
                    xl=self.params.xl,
                    xu=self.params.xu,
                    ftol=self.params.ftol,
                )
        except bsm.GSLException as exc:
            with self.failures.get_lock():
                self.failures.value += 1
            if self.failures.value >= self.Nfail:
                raise exc
            click.secho(f"{exc} raised", color="red")
            modes = np.array([np.nan] * 3)

        return [
            {"q": q, "theta": theta, "omega": m, "i": i}
            for (i, m) in enumerate(sorted(modes))
        ]


def data(fname, qs, thetas, params, hamiltonian=None, Nfail=int()):
    if hamiltonian is None:
        hamiltonian = params.cls == bsm.ModePolariton

    p = build_polariton(**params.as_args())

    qs *= p.state.delta / bsm.C

    pool = Pool()
    runner = Runner(p=p, hamiltonian=hamiltonian, params=params, Nfail=Nfail)
    with open(fname, "w") as f:
        writer = csv.DictWriter(f, ["q", "theta", "omega", "i"])
        writer.writeheader()
        for rows in tqdm.tqdm(
            pool.imap_unordered(runner, it.product(thetas, qs)),
            total=len(thetas) * len(qs),
        ):
            writer.writerows(rows)
