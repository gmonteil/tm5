#!/usr/bin/env python
from tkinter import OptionMenu

from flask import Flask, request, jsonify
import subprocess
from typing import Tuple
from numpy import base_repr
from loguru import logger
from pathlib import Path
import yaml
import os

app = Flask(__name__)


file_counter = (base_repr(i, base=36).lower() for i in range(360, 10_000_000))

def gencmd() -> Tuple[str, str]:
    run_id = next(file_counter)
    top_path = f'/data/avengers/fit_ic/results_flask-v2'
    while (Path(outpath := f'{top_path}/{run_id}').exists()):
        run_id = next(file_counter)
    datapath = f'/data/avengers/fit_ic/4server/current'
    fwd = f'{datapath}/forward.py'

    config = yaml.safe_load(request.form['conf'])
    emfile = config['emis']
    task=config['task']
    print(config)
    cmd = (
        '/home/guillaume/miniforge3/envs/tm5/bin/python '
        f'{fwd} --task {task} --output {outpath} '
        f'--emis {emfile} --data {datapath}'
    )
    logger.debug(cmd)
    return cmd, outpath


@app.route('/forward', methods=['POST'])
def run_forward():
    cmd, output_file = gencmd()
    logger.info(cmd)
    result = subprocess.run(cmd.split(), capture_output=False, check=True)
    return jsonify({'output': output_file})


if __name__ == "__main__":
    app.run(port=5018)
