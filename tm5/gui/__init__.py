# This could (should) be set through the yaml file or through an env variable. But keeping this here until we have a better mechanism
def get_hostname():
    import socket
    hostname = socket.gethostname()
    return hostname

# print(f"***{get_hostname().find('cosmos')}***")
import os
if get_hostname().find('cosmos')>=0:
    emission_dir = '/lunarc/nobackup/projects/ghg_inv/michael/TM5/input/ch4/emissions_fitic'
elif os.environ['TM5_HOST'] == 'gmlaptop':
    emission_dir = '/home/guillaume/iLab/TM5/emissions/CH4'
else:
    #-- on ICOS jupyter lab
    emission_dir = '/project/fit_ic/data/input/emissions/CH4'
    emission_dir = '/home/guillaume/iLab/TM5/emissions/CH4'


# Not sure what needs this, maybe we can delete it?
species_implemented = ['CO2', 'CH4']
