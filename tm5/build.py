#!/usr/bin/env python3
import os
from omegaconf import DictConfig
from pathlib import Path
import filecmp
import shutil
from typing import Dict, List
import tempfile
import subprocess
from omegaconf import OmegaConf
from tm5.run import run_tm5


def build_tm5(conf: DictConfig) -> Path:
    # Make build directory (and parents, if needed)
    Path(conf.build.directory).mkdir(exist_ok=True, parents=True)

    # Copy/Generate source files
    files = copy_files(conf.build)
    mfiles = gen_macros(conf.build)

    # Remove any source file that would no longer be needed
    cleanup(Path(conf.build.directory), list(files) + mfiles)

    # Create the makefile
    makefile = gen_makefile(conf.build)

    # Build TM5
    return make_tm5(makefile)


def copy_files(params : DictConfig) -> Dict[str, Path]:
    """
    Compile TM5
    - copy the source files in the build directory
    - create the makefile (with makedepf90)
    - build the executable
    - copy it to the run directory
    :return: path to the compiled executable
    """

    # Ensure that the build directory exists
    build_dir = Path(params.directory)
    build_dir.mkdir(exist_ok=True)

    # Make a list of all source files to be used:
    files = dict()
    for proj in params.paths.src:
        for file in Path(proj).iterdir():
            if file.is_file():
                files[file.name] = file

    # Copy the files to the build directory, if they differ from the files that are already in it:
    for filename, source in list(files.items()):
        dest = build_dir / filename

        # Remove the part after '__'? (default True)
        if params.paths.get("remove__", True):
            if '__' in filename :
                dest = build_dir / (filename.split('__')[0] + '.' + filename.split('.')[-1])
                files[dest.name] = filename

        if not dest.exists():
            shutil.copy(source, dest)
        elif not filecmp.cmp(source, dest):
            shutil.copy(source, dest)

    return files


def gen_macros(params: DictConfig) -> List[str]:
    """
    Generate files (*.inc files or sourcecode) that cannot be simply copied from source directories
    """

    # include files:
    fnames = []
    for fname, keys in params.macros.items():

        # Write the info to a temporary file
        tmpfile = ...
        fid, path = tempfile.mkstemp(text=True)
        with open(path, 'w') as fid:
            for key in keys :
                fid.write(f'#define {key}\n')

        # Copy the file to its final path only if it differs from the existing one:
        dest = Path(params.directory) / fname
        if not dest.exists():
            shutil.move(fid.name, dest)
        elif not filecmp.cmp(dest, fid.name):
            shutil.move(fid.name, dest)

        # Add the dest file to the list of files regardless
        fnames.append(dest.name)
    return fnames


def cleanup(build_dir: Path, files: List[str]):
    """
    Remove files from the build directory that :
    - do not exist in any of the source directory
    - were not created by gen_macros
    - were not created by an earlier compilation (*.mod and *.o and *.x files)
    :param params:
    :return:
    """
    for f in build_dir.iterdir():
        if f.name not in files and f.suffix not in ['.o', '.mod', '.x']:
            f.unlink()
        if f.name in ['dims_grid.o', 'dims_grid.mod']:
            f.unlink()


def gen_makefile(config: DictConfig) -> Path:
    """
    Create the Makefile_deps file and link to the correct makefile depending on make settings
    """

    # Choose the correct Makefile
    buildpath = Path(config.directory)
    makefile_source = buildpath / f'Makefile.{config.make.machine}.{config.make.compiler}.{config.make.options}'
    makefile_dest = buildpath / 'Makefile'
    shutil.copy(makefile_source, makefile_dest)

    # Append it with makedepf90
    # Temporarily move to the build directory for that
    curdir = os.getcwd()
    os.chdir(buildpath)
    command = ['makedepf90', '-o', 'tm5.x'] + list(Path().glob('*.[Ff]*'))
    with open(makefile_dest.name, 'a') as fid:
        subprocess.run(command, stdout=fid)
    os.chdir(curdir)

    return makefile_dest


def make_tm5(makefile : Path) -> Path:
    run_tm5(f'make -C {makefile.parent} {makefile.name} tm5.x'.split())
    return makefile.parent / 'tm5.x'


if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('--build', action='store_true', default=False, help='compile tm5')
    p.add_argument('config', nargs=1)
    args = p.parse_args()

    conf = OmegaConf.load(args.config)
    if args.build :
        build_tm5(conf)

