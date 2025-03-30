#!/usr/bin/env python


from flask import Flask, request, jsonify
from tm5.cli import workflows


# http://pancake.nebula:5000/forward -H "Content-Type: application/json" -d '{"command": "config_file"}'


app = Flask(__name__)


@app.route("/forward")
def run_command():
    args = request.json.get("command").split(' ')
    workflows.forward(args)