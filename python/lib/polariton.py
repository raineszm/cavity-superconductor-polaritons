import bardasis_schrieffer as bsm

from materials import NIOBIUM


def build_polariton(m, mu, Tc, vs, theta_v, T):
    sys = bsm.System(m, mu, Tc, vs, theta_v)
    state = bsm.State.solve(sys, T)
    bs = bsm.BS(state.delta/10, state)
    return bsm.Polariton(bs, bsm.Cavity(bs.root), bsm.Coupling(state))


def niobium(vs, theta_v=0., Trel=0.5):
    return build_polariton(vs=vs, T=Trel*NIOBIUM['Tc'], theta_v=theta_v, **NIOBIUM)
