from . import _bardasis_schrieffer as bsm

from .materials import NIOBIUM


def build_polariton(
    m, mu, Tc, vs, theta_s, T, root_rel, paraX, dipoleX, cls=bsm.Polariton
):
    sys = bsm.System(m, mu, Tc, vs, theta_s)
    state = bsm.State.solve(sys, T)
    bs = bsm.BS(0.3, state)
    return cls(
        coupling=bsm.Coupling(cav=bsm.Cavity(bs.root * root_rel), bs=bs),
        paraX=paraX,
        dipoleX=dipoleX,
    )


def niobium(
    vs, theta_s=0., Trel=0.5, root_rel=1., paraX=1., dipoleX=1., cls=bsm.Polariton
):
    return build_polariton(
        vs=vs,
        T=Trel * NIOBIUM["Tc"],
        theta_s=theta_s,
        root_rel=root_rel,
        paraX=paraX,
        dipoleX=dipoleX,
        cls=cls,
        **NIOBIUM
    )
