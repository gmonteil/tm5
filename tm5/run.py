#!/usr/bin/env python

import subprocess
from typing import List, Dict, Union
from loguru import logger
from pathlib import Path


def run_tm5(cmd: Union[str, List[str], List[Path]], output : str = None, settings : Dict = None):
    if isinstance(cmd, str):
        cmd = cmd.split()
    command = f'singularity run --cleanenv --bind {settings["tm5_path"]}:/tm5'
    command += f' {settings["container"]} '
    cmd = command.split() + cmd
    logger.info(' '.join([str(_) for _ in cmd]))
    try :
        return subprocess.run(cmd)
    except Exception as e:
        logger.exception(e)
        raise
