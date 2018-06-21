import click
import enum
from . import data
from .params import Params
from .ranges import Range
from pathlib import Path
import pendulum
import numpy as np
from .. import _bardasis_schrieffer as bsm


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
def main(qs, thetas, method, root_rel, dipole, para):
    m = Method[method]

    if m == Method.ACTION:
        cls = bsm.Polariton
    else:
        cls = bsm.ModePolariton

    params = Params(cls=cls, root_rel=root_rel, dipoleX=dipole, paraX=para)

    hamiltonian = m == Method.HAMILTONIAN

    fname = pendulum.now().format("MM-DD-YY_HH:MM:SS")
    fpath = Path(f"notebooks/data/{fname}_{method}.csv")

    data.data(str(fpath), qs, thetas, params, hamiltonian)

