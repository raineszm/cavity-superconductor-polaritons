"""
Material parameters for various known superconductors.

Units are converted to Gaussian atomic units
"""

from .units import ureg, Q_

# Atomic units in meV
"""Material parameters for Niobium """
NIOBIUM = dict(  # meV
    Tc=Q_(9.5, "kelvin").to("meV", "boltzmann").magnitude,
    mu=Q_(6.18e4, "kelvin").to("meV", "boltzmann").magnitude,
    m=Q_(1.6, "m_e").magnitude,  # bare electron mass
)
