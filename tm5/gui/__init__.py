from omegaconf import OmegaConf
from pathlib import Path
host = OmegaConf.load(Path.home() / ".conf/fitic/gui.conf")