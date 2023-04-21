#!/usr/bin/env python

from dataclasses import dataclass
from pint import Unit, Quantity
from lumia.units import units_registry as units
from lumia.units import xmc, xmo, xm14c


@dataclass
class Specie:
    unit_emis   : Unit
    unit_mix    : Unit
    molar_mass  : Quantity

    def __post_init__(self):
        self.unit_mix = units(self.unit_mix).units
        self.unit_emis = units(self.unit_emis).units

    def __hash__(self):
        # This is a workaround for python>=3.11, where dataclasses only accept hashable objects as defaults (but don't actually check the hash)
        ...


species = {
    'CO2' : Specie(unit_emis = 'PgCO2', unit_mix = 'ppm', molar_mass= Quantity(f'{xmc + 2 * xmo} g/mol')),
    '14CO2' : Specie(unit_emis = 'PgCO2', unit_mix = 'ppm', molar_mass= Quantity(f'{xm14c + 2 * xmo} g/mol'))
}
