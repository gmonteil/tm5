from omegaconf import OmegaConf
from pathlib import Path
import os

try:
    # The following should be set by the conda env. If not, use (once) 1conda env vars set TM5_HOST=my_host'
    # and there should be a corresponding "my_host" section in gui.conf
    
    hostname = os.environ['TM5_HOST']
    host = OmegaConf.load(Path.home() / ".config/fitic/gui.conf")[hostname]
except FileNotFoundError:
    host = None
