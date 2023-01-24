#!/usr/bin/env python

"""
Main run scenarios
"""

from omegaconf import OmegaConf, DictConfig


def forward(config: DictConfig):
    """
    Run a TM5 forward simulation
    The config file is assumed to have been filled in already, and contains all the required information
    """
    # 1. Build the observations file (if any), and other auxiliary files (stations list, etc.)

    # 2. Build the emissions file

    # 3. Copy the initial condition

    # 4. Mount the container

    # 5. Create the rc-file

    # 6. Run

    # 7. Post-process


if __name__ == '__main__':
    pass