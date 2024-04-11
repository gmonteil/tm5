#!/usr/bin/env python

from dataclasses import dataclass
from pint import Unit
from tm5.units import units_registry as units


@dataclass
class Specie:
    unit_emis   : Unit
    unit_mix    : Unit

    def __post_init__(self):
        self.unit_mix = units(self.unit_mix).units

    def __hash__(self):
        # This is a workaround for python>=3.11, where dataclasses only accept hashable objects as defaults (but don't actually check the hash)
        ...


species = {
    'CO2' : Specie(unit_emis = 'PgCO2', unit_mix = 'ppm'),
    '14CO2' : Specie(unit_emis = 'PgCO2', unit_mix = 'ppm'),
    'CH4' : Specie(unit_emis = 'PgCH4', unit_mix = 'ppb')
}
