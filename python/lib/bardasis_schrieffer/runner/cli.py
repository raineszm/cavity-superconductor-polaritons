import enum
import json
from pathlib import Path

import click
import numpy as np
import pendulum

from . import data
from .. import _bardasis_schrieffer as bsm
from .notify import push_notification
from .params import Params
from .ranges import Range


class Method(enum.Enum):
    ACTION = enum.auto()
    MODEACTION = enum.auto()
    HAMILTONIAN = enum.auto()


@click.command()
@click.option("--qs", type=Range(), default="0:0.4:100")
@click.option("--thetas", type=Range(), default=np.array([0.]))
@click.option(
    "--method", type=click.Choice(list(Method.__members__)), default="HAMILTONIAN"
)
@click.option("-r", "--root_rel", type=float, default=1.)
@click.option("--dipole", type=float, default=1.)
@click.option("--para", type=float, default=1.)
@click.option("--notify/--no-notify", default=False)
def main(qs, thetas, method, root_rel, dipole, para, notify):
    m = Method[method]

    if m == Method.ACTION:
        cls = bsm.Polariton
    else:
        cls = bsm.ModePolariton

    params = Params(cls=cls, root_rel=root_rel, dipoleX=dipole, paraX=para)

    hamiltonian = m == Method.HAMILTONIAN

    fname = pendulum.now().format("MM-DD-YY_HH:MM:SS")
    fpath_root = Path(f"notebooks/data/{fname}_{method}")

    meta_args = dict(method=method, params=params.as_args())
    del meta_args["params"]["cls"]

    with open(fpath_root.with_suffix(".json"), "w") as json_file:
        json.dump(meta_args, json_file)

    data.data(fpath_root.with_suffix(".csv"), qs, thetas, params, hamiltonian)
    if notify:
        push_notification(fname, len(qs) * len(thetas))
