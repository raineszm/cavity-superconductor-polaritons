from dataclasses import dataclass, asdict, field
from typing import Mapping, Type
from ..materials import *
from .. import _bardasis_schrieffer as bsm


@dataclass(frozen=True)
class Params:
    material: Mapping[str, float] = field(default_factory=lambda: PNICTIDE)
    Trel: float = 0.4
    vrel: float = 0.6
    theta_s: float = 0.
    dipoleX: float = 1
    paraX: float = 1
    root_rel: float = 1
    cls: Type[bsm.Polariton] = bsm.ModePolariton

    def as_args(self):
        pdict = asdict(self)
        vrel = pdict.pop("vrel")
        pdict["vs"] = vrel * self.vc()
        return pdict

    def vc(self):
        sys0 = bsm.System(vs=0., **self.material)
        state0 = bsm.State.solve(sys0, self.Trel * self.material["Tc"])
        return (
            state0.delta / sys0.kf
        )  # The critical superfluid velocity via Landau's argument

