#!/usr/bin/env python
from pint import UnitRegistry

units_registry = UnitRegistry()

xmc = 12.01115
xmo = 15.9994

units_registry.define(f'gC = {1 / xmc} mol')
units_registry.define(f'gCO2 = {1 / (xmc + 2*xmo)} mol')

units_registry.define('ppm = umol/mol')
units_registry.define('ppb = nmol/mol')