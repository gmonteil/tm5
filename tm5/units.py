#!/usr/bin/env python

from pint import UnitRegistry

units_registry = UnitRegistry()

xmh   =  1.00784
xmc   = 12.01115
xm14c = 14.003241989
xmo   = 15.9994

units_registry.define(f'gH = {1 / xmh} mol')
units_registry.define(f'gC = {1 / xmc} mol')
units_registry.define(f'gCO2 = {1 / (xmc + 2*xmo)} mol')
units_registry.define(f'g14CO2 = {1 / (xm14c + 2*xmo)} mol')
units_registry.define(f'gCH4 = {1 / (xmc + 4*xmh)} mol')

units_registry.define('ppm = umol/mol')
units_registry.define('ppb = nmol/mol')

if __name__=="__main__":
    ureg = units_registry
    # print(f"{ureg.metric_ton/1000}")
    # print(f"{ureg}")
    # print(f"{ureg.tonne}")
    mass_lst = ['kg', 'g'] #-- does not work with tonne ['t','tonne']
    species_lst = ['CO2','CH4', 'H', 'C']
    for species in species_lst:
        for m in mass_lst:
            value = ureg.Quantity('mol').to(f'{m}{species}').m
            print(f"1 mol of {species:>5s}: {value:.6E} [{m}]")
    for species in species_lst:
        value = ureg.Quantity('kg{species}').to('mol')
        print(f"1kg of {species:>5s}: {value:.6E} [mol]")
    #
    # species = 'H'
    # print(ureg.Quantity('mol').to('tH').m)
    # print(ureg.Quantity('mol').to('tonH').m)
    # print(ureg.Quantity('mol').to('tonneH').m)
    for m in ['t','kt','Mt','Gt']:
        value = ureg.Quantity('mol').to(f'{m}{species}').m
        print(f"1 mol of {species:>5s}: {value:.6E} [{m}]")
