from units import ureg, Q_

# Atomic units in meV
NIOBIUM = dict(  # meV
    Tc=Q_(9.5, 'kelvin').to('meV', 'boltzmann').magnitude,
    mu=Q_(6.18e4, 'kelvin').to('meV', 'boltzmann').magnitude,
    m=1.6  # bare electron mass
)
