from . import _bardasis_schrieffer as bsm

from .materials import NIOBIUM


def _build_polariton(
    m, mu, Tc, vs, theta_s, Trel, root_rel, paraX, dipoleX, cls=bsm.Polariton
):
    sys = bsm.System(m, mu, Tc, vs, theta_s)
    state = bsm.State.solve(sys, Trel * Tc)
    bs = bsm.BS(0.3, state)
    return cls(
        coupling=bsm.Coupling(cav=bsm.Cavity(bs.root * root_rel), bs=bs),
        paraX=paraX,
        dipoleX=dipoleX,
    )


def build_polariton(material=NIOBIUM, **kwargs):
    return _build_polariton(**material, **kwargs)
