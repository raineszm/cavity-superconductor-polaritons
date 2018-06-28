import csv
import itertools as it
from multiprocessing import Pool, Value, Manager
import logbook
import logbook.queues as lq
from pathlib import Path

import click
import numpy as np
import tqdm

from .. import _bardasis_schrieffer as bsm
from ..materials import NIOBIUM
from ..polariton import build_polariton


class Runner:
    def __init__(self, p, hamiltonian, params, Nfail, failures, queue):
        self.p = p
        self.hamiltonian = hamiltonian
        self.params = params
        self.Nfail = Nfail
        self.failures = failures
        self.handler = lq.MultiProcessingHandler(queue)

    def __call__(self, args):
        (theta, q) = args
        try:
            if self.hamiltonian:
                modes = self.p.bands(q, theta)
            else:
                modes = self.p.find_modes(
                    q,
                    theta,
                    xl=self.params.xl * self.p.bs.root,
                    xu=self.params.xu * self.p.bs.root,
                    ftol=self.params.ftol,
                )
        except bsm.GSLException as exc:
            self.failures.value += 1
            if self.failures.value >= self.Nfail:
                raise exc
            click.secho(f"{exc} raised", color="red")
            with self.handler.applicationbound():
                log = logbook.Logger('Runner')
                log.error(exc)
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

    manager = Manager()
    subscriber = lq.MultiProcessingSubscriber(manager.Queue(-1))
    logfile = logbook.FileHandler(fname.with_suffix(".log"), delay=True)

    runner = Runner(
        p=p,
        hamiltonian=hamiltonian,
        params=params,
        Nfail=Nfail,
        failures=manager.Value("i", 0),
        queue=subscriber.queue,
    )
    controller = subscriber.dispatch_in_background(logfile)

    with open(fname, "w") as f:
        writer = csv.DictWriter(f, ["q", "theta", "omega", "i"])
        writer.writeheader()

        pool = Pool()
        for rows in tqdm.tqdm(
            pool.imap_unordered(runner, it.product(thetas, qs)),
            total=len(thetas) * len(qs),
        ):
            writer.writerows(rows)
    controller.stop()
