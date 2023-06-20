#!/usr/bin/env python

import subprocess
from typing import List, Dict, Union
from loguru import logger
from pathlib import Path


def run_tm5(cmd: str | List[str] | List[Path], settings : Dict = None, stdout = None):
    """
    Run a command, possibly inside a singularity container (preconfigured for TM5).
    Arguments :
        * cmd : the command to run (can be a string or a list of strings)
        * settings : should be a dictionary (or an omegaconf.DictConfig) with the following keys:
            - container: path towards the singularity container (*.simg file). If the "container" key is not provided, then the command "cmd"
             is just run as a subprocess.
            - tm5_path: points toward the path to be mounted under "/tm5"
            - bind (optional): contains a list of additional paths of the host that should be mounted inside the container.

    Example :
    The following line will run the "ipython" command in a container located at "/home/user/TM5/images/tm5.sif", and with the "/home/user/TM5", "/data" and "/scratch" of the host mounted respectively under "/tm5", "/data" and "/scratch" in the container.

    _ = run_tm5('ipython', settings = {
                    'tm5_path': '/home/user/TM5',
                    'container': '/home/user/TM5/images/tm5.sif',
                    'bind': ['/data', '/scratch']})

    By contrast, the command "run_tm5('ipython', settings={'tm5_path':'/home/user/TM5'})" will just spawn a ipython shell in a subprocess, since no "container" key is provided in "settings". The other keys in "settings" will be ignored.
    """

    # Ensure that the command is given as a list of strings
    if isinstance(cmd, str):
        cmd = cmd.split()

    # if 'container' in settings :
    if settings is not None:
        # Prepare the singularity-related part of the command
        command = f'singularity run --cleanenv --bind {settings["tm5_path"]}:/tm5'
        for mountpath in settings.get('bind', []):
            command += f' --bind {mountpath}'
        command += f' {settings["container"]} '

        # Append the actual command, log and run:
        cmd = command.split() + cmd

    logger.info(' '.join([str(_) for _ in cmd]))
    try :
        _ = subprocess.run(cmd, stdout=stdout)
        if _.returncode != 0 :
            raise RuntimeError(_)
        return _
    except Exception as e:
        logger.exception(e)
        raise
