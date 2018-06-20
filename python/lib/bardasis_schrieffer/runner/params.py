from dataclasses import dataclass, asdict, field
from typing import Dict, Type
from ..materials import NIOBIUM
from .. import _bardasis_schrieffer as bsm


@dataclass
class Params:
    material: Dict[str, float] = field(default_factory=lambda: NIOBIUM)
    Trel: float = 0.4
    vrel: float = 0.6
    theta_s: float = 0.
    dipoleX: float = 1
    paraX: float = 1
    root_rel: float = 1
    cls: Type[bsm.Polariton] = bsm.ModePolariton

    def as_args(self, vc):
        pdict = asdict(self)
        vrel = pdict.pop("vrel")
        pdict["vs"] = vrel * vc
        return pdict
