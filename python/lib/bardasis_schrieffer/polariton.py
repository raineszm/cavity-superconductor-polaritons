from . import _bardasis_schrieffer as bsm

from .materials import NIOBIUM


def _build_polariton(
    m, mu, Tc, vs, theta_s, Trel, root_rel, paraX, dipoleX, bsmass, cls=bsm.Polariton
):
    sys = bsm.System(m, mu, Tc, vs, theta_s)
    state = bsm.State.solve(sys, Trel * Tc)
    bs = bsm.BS(bsmass, state)
    sys0 = bsm.System(m, mu, Tc, vs, 0)
    state0 = bsm.State.solve(sys0, Trel * Tc)
    bs0 = bsm.BS(bsmass, state0)
    return cls(
        coupling=bsm.Coupling(cav=bsm.Cavity(bs0.root * root_rel), bs=bs),
        paraX=paraX,
        dipoleX=dipoleX,
    )


def build_polariton(material=NIOBIUM, **kwargs):
    return _build_polariton(**material, **kwargs)
