import bardasis_schrieffer as bsm

NIOBIUM = dict(
    Tc=3e-5,  # Hartrees
    mu=0.1955,  # Hartress
    m=1.6  # bare electron
)


def build_polariton(m, mu, Tc, vs, theta_v, T):
    sys = bsm.System(m, mu, vs, theta_v)
    state = bsm.State.solve(sys, T)
    bs = bsm.BS(state.delta, state)
    return bsm.Polariton(bs, bsm.Cavity(bs.root), bsm.Coupling(state))


def niobium(vs):
    return build_polariton(vs=vs, T=NIOBIUM['Tc']/2, theta_v=0., **NIOBIUM)
