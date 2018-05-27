import bardasis_schrieffer as bsm

from materials import NIOBIUM


def build_polariton(m, mu, Tc, vs, theta_v, T, big=1):
    sys = bsm.System(m, mu, Tc, vs, theta_v)
    state = bsm.State.solve(sys, T)
    bs = bsm.BS(state.delta/10, state)
    return bsm.Polariton(bs, bsm.Cavity(bs.root*0.8), bsm.Coupling(state), big=big)


def niobium(vs, theta_v=0., Trel=0.5, big=1.):
    return build_polariton(vs=vs, T=Trel*NIOBIUM['Tc'], theta_v=theta_v, big=big, **NIOBIUM)
