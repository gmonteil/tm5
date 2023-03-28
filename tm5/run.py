#!/usr/bin/env python

import subprocess
from typing import List


def run_tm5(cmd: List[str], output : str = None):
    options = []
    if output is not None :
        options = ['--output', output]
    cmd = ['tm5', '--dev'] + options + cmd
    print(cmd)
    return subprocess.run(cmd)