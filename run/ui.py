"""
Very simple interface defining:
- a settings window/widget ==> simple yaml file editor
- a button to build/run the code
- a few toogles to enable run options (e.g. build or not)
- a terminal to display the code running
- a tab showing some diagnostics plots
"""

import panel as pn
from types import SimpleNamespace
from loguru import logger
import sys, os
import tm5
import tempfile


pn.extension()
pn.extension('codeeditor')
pn.extension('terminal')


#=============================================================
# User settings (might be replaced by an ArgumentParser)
args = SimpleNamespace(
    config='forward.yaml',
    host = os.environ['TM5_HOST']
)
logger.info(args)


#=============================================================
# Widgets container
widgets = {}


#=============================================================
# Widget creation functions
def create_editor(fname: str, sizing_mode: str = 'stretch_both', language: str='yaml'):
    with open(fname, 'r') as fid:
        lines = ''.join(fid.readlines())
    return pn.widgets.CodeEditor(
        value=lines, 
        sizing_mode=sizing_mode, 
        language=language,
        )


#=============================================================
# Action functions
#def save_yaml(event: bool):
#    if not event : return
#    term = widgets['terminal']
#    with open(args.config, 'w') as fid :
#        fid.writelines(widgets['editor'].value)
#    term.clear()
#    term.subprocess.run('cat', args.config)


def save_yaml(filename: str | None = None) -> str:
    """
    Write the content of the text editor to a yaml file (either provided as an argument, or, by default, to a temporary file), and return the path to that file
    """

    if filename is None:
        fid, filename = tempfile.mkstemp()
        os.close(fid)
    with open(filename, 'w') as fid:
        fid.writelines(widgets['editor'].value)
    return filename


def run_TM5(event: bool):
    if not event: return
    widgets['terminal'].subprocess.run('python', 'forward.py', '-b', '-m', args.host, save_yaml())


#=============================================================
# Widget definitions 
widgets['editor'] = create_editor(args.config)
widgets['run_button'] = pn.widgets.Button(name='Run TM5', button_type='primary')
widgets['terminal'] = pn.widgets.Terminal(options={"cursorBlink": True}, height=600, sizing_mode='stretch_width')


#=============================================================
# Bind actions to widgets
#pn.bind(save_yaml, widgets['run_button'], watch=True)
widgets['run_button'].on_click(run_TM5)


#=============================================================
# Global layout
template = pn.template.MaterialTemplate(title='FIT-IC')
template.main.append(
    pn.Tabs(('yaml editor', pn.Row(
        widgets['editor'],
        pn.Column(
            widgets['run_button'],
            widgets['terminal']
            )
        )
    ))
)
template.servable()