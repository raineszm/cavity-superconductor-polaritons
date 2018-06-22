"""
Material parameters for various known superconductors.

Units are converted to Gaussian atomic units
"""

from .units import ureg, Q_

SCALAR = Q_(1, "hartree").to("meV").magnitude

# Atomic units in meV
"""Material parameters for Niobium """
NIOBIUM = dict(  # meV
    Tc=Q_(9.5, "kelvin").to("meV", "boltzmann").magnitude,
    mu=Q_(6.18e4, "kelvin").to("meV", "boltzmann").magnitude,
    m=Q_(1.6, "m_e").magnitude * SCALAR,  # bare electron mass * meV/Hartree
)

"""Material parameters for a model of pnictides"""
PNICTIDE = dict(  # meV
    Tc=Q_(35, "kelvin").to("meV", "boltzmann").magnitude,
    mu=Q_(100, "meV").magnitude,
    m=Q_(0.7, "m_e").magnitude * SCALAR,  # bare electron mass * meV/Hartree
)
