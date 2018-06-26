"""
Material parameters for various known superconductors.

Units are converted to Gaussian atomic units
"""

from .units import ureg, Q_

"""Material parameters for Niobium """
NIOBIUM = dict(  # eV
    Tc=Q_(9.5, "kelvin").to("eV", "boltzmann").magnitude,
    mu=Q_(6.18e4, "kelvin").to("eV", "boltzmann").magnitude,
    m=Q_(1.6, "m_e*c^2").to("eV").magnitude,
)

"""Material parameters for a model of pnictides"""
PNICTIDE = dict(  # eV
    Tc=Q_(35, "kelvin").to("eV", "boltzmann").magnitude,
    mu=Q_(100, "meV").to("eV").magnitude,
    m=Q_(0.7, "m_e*c^2").to("eV").magnitude,
)
