from .tracers import CH4TracerSettings, CO2TracerSettings
from .emissions import EmissionSettings
from .runsettings import RunSettings
from .reactions import ReactionSettings
import panel as pn


# This could (should) be set through the yaml file or through an env variable. But keeping this here until we have a better mechanism
def get_hostname():
    import socket
    hostname = socket.gethostname()
    return hostname

# print(f"***{get_hostname().find('cosmos')}***")
if get_hostname().find('cosmos')>=0:
    emission_dir = '/lunarc/nobackup/projects/ghg_inv/michael/TM5/input/ch4/emissions_fitic'
else:
    #-- on ICOS jupyter lab
    emission_dir = '/project/fit_ic/data/input/emissions/CH4'
    
    
# Not sure what needs this, maybe we can delete it?
species_implemented = ['CO2', 'CH4']