import bardasis_schrieffer as bsm

from materials import NIOBIUM


def build_polariton(m, mu, Tc, vs, theta_v, T, root_rel, big):
    sys = bsm.System(m, mu, Tc, vs, theta_v)
    state = bsm.State.solve(sys, T)
    bs = bsm.BS(0.3, state)
    return bsm.Polariton(bs, bsm.Cavity(bs.root*root_rel), bsm.Coupling(state), big)


def niobium(vs, theta_v=0., Trel=0.5, root_rel=1., big=1.):
    return build_polariton(vs=vs, T=Trel*NIOBIUM['Tc'], theta_v=theta_v, root_rel=root_rel, big=big, **NIOBIUM)
