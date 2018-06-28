import enum
import json
import sys
from pathlib import Path

import click
import numpy as np
import pendulum

from .. import __version__
from . import data
from .. import _bardasis_schrieffer as bsm
from .notify import notify_failure, notify_success
from .params import Params
from .ranges import Range


class Method(enum.Enum):
    ACTION = enum.auto()
    MODEACTION = enum.auto()
    HAMILTONIAN = enum.auto()


def show_param(name, value):
    if hasattr(value, "__len__"):
        return "{}: [{}, {}] x {}".format(name, min(value), max(value), len(value))
    else:
        return "{}: {}".format(name, value)


def get_confirmation(**params):
    click.secho("Please confirm parameters:", fg="green")
    for (k, v) in params.items():
        click.echo(show_param(k, v))
    click.confirm("Is this correct?", abort=True)


@click.command()
@click.version_option(__version__)
@click.option("--qs", type=Range(), default="0:0.4:100")
@click.option("--thetas", type=Range(), default=np.array([0.]))
@click.option(
    "--method", type=click.Choice(list(Method.__members__)), default="HAMILTONIAN"
)
@click.option("--notify/--no-notify", default=False)
@click.option("--confirm/--no-confirm", default=True)
@click.option("--data-dir", type=click.Path(exists=True), default=None)
@click.option("--Nfail", type=int, default=int(1e8))
@click.option("-r", "--root_rel", type=float, default=Params().root_rel)
@click.option("--dipole", type=float, default=Params().dipole)
@click.option("--para", type=float, default=Params().para)
@click.option("-v", "--vrel", type=float, default=Params().vrel)
@click.option("--xl", type=float, default=Params().xl)
@click.option("--xu", type=float, default=Params().xu)
@click.option("--ftol", type=float, default=Params().ftol)
def main(qs, thetas, method, notify, confirm, data_dir, nfail, **kwargs):
    m = Method[method]

    if m == Method.ACTION:
        cls = bsm.Polariton
    else:
        cls = bsm.ModePolariton

    params = Params(cls=cls, **kwargs)

    hamiltonian = m == Method.HAMILTONIAN

    if hamiltonian:
        for k in Params.REMOVE:
            if k in kwargs:
                click.secho(f"{k} will be ignored with method {m}", color="red")

    if confirm:
        get_confirmation(hamiltonian=hamiltonian, **params.asdict())

    fname = pendulum.now().format("MM-DD-YY_HH:MM:SS")
    fpath_root = Path(f"notebooks/data/{fname}_{method}")

    meta_args = dict(
        version=__version__,
        method=method,
        params=params.as_args(),
        extra=dict(xl=params.xl, xu=params.xu, ftol=params.ftol),
    )
    del meta_args["params"]["cls"]

    jsonpath = fpath_root.with_suffix(".json")
    with open(jsonpath, "w") as json_file:
        json.dump(meta_args, json_file)

    datapath = fpath_root.with_suffix(".csv")
    try:
        data.data(
            fname=datapath,
            qs=qs,
            thetas=thetas,
            params=params,
            hamiltonian=hamiltonian,
            Nfail=nfail,
        )
    except bsm.GSLException as exc:
        if notify:
            notify_failure(datapath, nfail, exc)

        sys.exit(-1)

    if notify:
        notify_success(datapath, len(qs) * len(thetas))

    if data_dir:
        import shutil

        shutil.copy(datapath, data_dir)
        shutil.copy(jsonpath, data_dir)
