#!/usr/bin/env python
from tkinter import OptionMenu

import sys
from flask import Flask, request, jsonify
import subprocess
from typing import Tuple
from numpy import base_repr
from loguru import logger
from pathlib import Path
import yaml
import f90nml
import os
import shutil
from argparse import ArgumentParser

app = Flask(__name__)

#logger.add('/srv/data/AVENGERS/fit_ic/gui/logs/fitic-flask-server.log', level='DEBUG')

file_counter = (base_repr(i, base=36).lower() for i in range(360, 10_000_000))

# parser = ArgumentParser()
# parser.add_argument('datapath',
#                     type=Path,
#                     default='/data/avengers/fit_ic/4server/current',
#                     help='Path where inversion inputs are found')
# parser.add_argument('resultdir',
#                     type=Path,
#                     default='/data/avengers/fit_ic/results_flask-v3',
#                     help='top-level result directory, the individual user output directories will be generated here.')

# args = parser.parse_args(sys.argv[1:])

# datapath = args.datapath
# resultdir = args.resultdir

#
#-- ATTENTION: path settings depend on platform!!!
#
resultdir = f'/data/avengers/fit_ic/results_flask-v3'
python = '/data/avengers/python/fitic/bin/python'
datapath = f'/data/avengers/fit_ic/4server/current'

def get_host_port():
    #MVO-ATTENTION: does not work properly on pancake...
    host = os.getenv("FLASK_RUN_HOST", "127.0.0.1") # 127.0.0.1 is default value when FLASK_RUN_HOST is not set
    port = os.getenv("FLASK_RUN_PORT", "5000") # 5000 is default value when FLASK_RUN_PORT is not set
    #print(f"App is running at http://{host}:{port}") #Btw it should be "App is scheduled to run at..."
    return host,port

def gencmd() -> Tuple[str, str]:
    # _, port = get_host_port()
   
    run_id = next(file_counter)
    while (Path(outpath := f'{resultdir}/{run_id}').exists()):
        run_id = next(file_counter)

    fwd = f'{datapath}/forward.py'
    config = yaml.safe_load(request.form['conf'])
    emfile = config['emis']
    task = config['task']
    logger.debug(f"config ***{config}***")
    # print(config)
    cmd = (
        f'{python} '
        f'{fwd} --task {task} --output {outpath} '
        f'--emis {emfile} --data {datapath}'
    )
    #
    #--
    #
    outpath = Path(outpath)
    outpath.mkdir(parents=True, exist_ok=True)
    if 'namelist' in config:
        nmlfile = outpath / 'fitic.nml'
        nml_dict = { 'fitic.nml' : config['namelist'] }
        nml = f90nml.Namelist(nml_dict)
        nml.write(nmlfile, force=True)

    logger.debug(f"command ==>{cmd}<==")

    return cmd, str(outpath)


@app.route('/forward', methods=['POST'])
def run_forward():
    cmd, output_file = gencmd()
    logger.info(cmd)
    result = subprocess.run(cmd.split(), capture_output=False, check=True)
    return jsonify({'output': output_file})



#===============================================================================
#
# OLD ROUTINES CURRENTLY NOT USED
#
#===============================================================================

@app.route('/submit')
def submit():
    cmd, output_file = 'tsp ' + gencmd()

    # Submit the run with tsp
    jobid = subprocess.run(cmd.split(), capture_output=True, text=True, check=True).stdout.strip()

    # Status will be one of "queued", "finished" or "running"
    status = subprocess.run(['tsp', '-s', jobid], capture_output=True, text=True, check=True).stdout.strip()

    outfile = ''
    if status != 'queued':
        # Name of the output file
        outfile = subprocess.run(['tsp', '-o', jobid], capture_output=True, text=True, check=True).stdout.strip()

    return jsonify({"jobid": jobid, 'status':status, 'outfile':outfile})


@app.route('/status/<jobid>')
def check_status(jobid):
    status = subprocess.run(['tsp', '-s', jobid], capture_output=True, text=True, check=True).stdout.strip()
    info = subprocess.run(['tsp', '-i', jobid], capture_output=True, text=True, check=True).stdout.strip()
    outfile = subprocess.run(['tsp', '-o', jobid], capture_output=True, text=True, check=True).stdout.strip()
    with open(outfile) as fid:
        outlines = fid.readlines()
    return jsonify({'status': status, 'info': info, 'outfile':outfile, 'stdout':outlines})


if __name__ == "__main__":
    app.run(port=5018)
