#!/usr/bin/env python

import subprocess
from typing import List, Dict, Union
from loguru import logger
from pathlib import Path


def run_tm5(cmd: Union[str, List[str], List[Path]], settings : Dict = None):
    """
    Run a command inside a singularity container (preconfigured for TM5).
    Arguments :
        * cmd : the command to run (can be a string or a list of strings)
        * settings : should be a dictionary (or an omegaconf.DictConfig) with the following keys:
            - tm5_path: points toward the path to be mounted under "/tm5"
            - container: path towards the singularity container (*.simg file)
            - bind (optional): contains a list of additional paths of the host that should be mounted inside the container.

    Example :
    The following line will run the "ipython" command in a container located at "/home/user/TM5/images/tm5.sif", and with the "/home/user/TM5", "/data" and "/scratch" of the host mounted respectively under "/tm5", "/data" and "/scratch" in the container.

    _ = run_tm5('ipython', settings = {
                    'tm5_path': '/home/user/TM5',
                    'container': '/home/user/TM5/images/tm5.sif',
                    'bind': ['/data', '/scratch']})
    """

    # Ensure that the command is given as a list of strings
    if isinstance(cmd, str):
        cmd = cmd.split()

    # Prepare the singularity-related part of the command
    command = f'singularity run --cleanenv --bind {settings["tm5_path"]}:/tm5'
    for mountpath in settings.get('bind', []):
        command += f' --bind {mountpath}'
    command += f' {settings["container"]} '

    # Append the actual command, log and run:
    cmd = command.split() + cmd
    logger.info(' '.join([str(_) for _ in cmd]))
    try :
        return subprocess.run(cmd)
    except Exception as e:
        logger.exception(e)
        raise
