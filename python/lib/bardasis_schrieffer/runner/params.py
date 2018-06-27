from dataclasses import dataclass, asdict, field
from typing import Mapping, Type, ClassVar, List
from ..materials import *
from .. import _bardasis_schrieffer as bsm


@dataclass(frozen=True)
class Params:
    material: Mapping[str, float] = field(default_factory=lambda: PNICTIDE)
    Trel: float = 0.4
    vrel: float = 0.6
    theta_s: float = 0.
    dipole: float = 1
    para: float = 1
    root_rel: float = 1
    bsmass: float = 0.1
    cls: Type[bsm.Polariton] = bsm.ModePolariton
    xl: float = 0.8
    xu: float = 1.2
    ftol: float = 1e-2

    RENAME: ClassVar[Mapping[str, str]] = {"dipole": "dipoleX", "para": "paraX"}
    REMOVE: ClassVar[List[str]] = ["xl", "xu", "ftol"]

    def as_args(self):
        pdict = asdict(self)

        for (k, v) in self.RENAME.items():
            pdict[v] = pdict.pop(k)

        vrel = pdict.pop("vrel")
        pdict["vs"] = vrel * self.vc()

        for k in self.REMOVE:
            del pdict[k]

        return pdict

    def vc(self):
        sys0 = bsm.System(vs=0., **self.material)
        state0 = bsm.State.solve(sys0, self.Trel * self.material["Tc"])
        return (
            state0.delta / sys0.kf
        )  # The critical superfluid velocity via Landau's argument

