'''
Unit conversion.

We use Gaussian atomic units for these calculations
'''
import pint
ureg = pint.UnitRegistry(system='cgs')
Q_ = ureg.Quantity
