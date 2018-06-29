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
    def __init__(self, hamiltonian, params, Nfail, failures, queue):
        self._p = build_polariton(**params.as_args())
        self.hamiltonian = hamiltonian
        self.params = params
        self.Nfail = Nfail
        self.failures = failures
        self.handler = lq.MultiProcessingHandler(queue)

    def p(self, theta_s):
        if theta_s != self._p.state.sys.theta_s:
            pdict = self.params.as_args()
            pdict["theta_s"] = theta_s
            self._p = build_polariton(**pdict)
        return self._p

    def __call__(self, args):
        (theta_s, theta, q) = args
        p = self.p(theta_s * np.pi)
        try:
            if self.hamiltonian:
                modes = p.bands(q, theta * np.pi)
            else:
                modes = p.find_modes(
                    q,
                    theta * np.pi,
                    xl=self.params.xl * p.bs.root,
                    xu=self.params.xu * p.bs.root,
                    ftol=self.params.ftol,
                )
        except bsm.GSLException as exc:
            self.failures.value += 1
            if self.failures.value >= self.Nfail:
                raise exc
            click.secho(f"{exc} raised", color="red")
            with self.handler.applicationbound():
                log = logbook.Logger("Runner")
                log.error(exc)
            modes = np.array([np.nan] * 3)

        return [
            {"theta_s": theta_s, "q": q, "theta": theta, "omega": m, "i": i}
            for (i, m) in enumerate(sorted(modes))
        ]


def data(fname, q, theta, theta_s, params, hamiltonian=None, Nfail=int()):
    if hamiltonian is None:
        hamiltonian = params.cls == bsm.ModePolariton

    p = build_polariton(**params.as_args())

    q *= p.state.delta / bsm.C

    manager = Manager()
    subscriber = lq.MultiProcessingSubscriber(manager.Queue(-1))
    logfile = logbook.FileHandler(fname.with_suffix(".log"), delay=True)

    runner = Runner(
        hamiltonian=hamiltonian,
        params=params,
        Nfail=Nfail,
        failures=manager.Value("i", 0),
        queue=subscriber.queue,
    )
    controller = subscriber.dispatch_in_background(logfile)

    with open(fname, "w") as f:
        writer = csv.DictWriter(f, ["theta_s", "q", "theta", "omega", "i"])
        writer.writeheader()

        pool = Pool()
        for rows in tqdm.tqdm(
            pool.imap_unordered(runner, it.product(theta_s, theta, q)),
            total=len(theta_s) * len(theta) * len(q),
        ):
            writer.writerows(rows)
    controller.stop()
