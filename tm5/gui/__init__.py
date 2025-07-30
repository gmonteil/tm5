from omegaconf import OmegaConf
from pathlib import Path
try:
    host = OmegaConf.load(Path.home() / ".config/fitic/gui.conf")
except FileNotFoundError:
    host = None
