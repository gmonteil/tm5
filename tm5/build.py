#!/usr/bin/env python3
import os
from omegaconf import DictConfig
from pathlib import Path
import filecmp
import shutil
from typing import Dict, List
import tempfile
from omegaconf import OmegaConf
from tm5.run import run_tm5
from loguru import logger
import subprocess
from tm5 import debug
from tm5.system import runcmd


@debug.trace_call
def build_tm5(dconf: DictConfig, clean : bool = False) -> Path:
    """
    Compile TM5, based on the settings passed through the "conf" argument (i.e. the main yaml file). The yaml file should have a "build" section with the following structure: 

    build :
        directory : path where the source code will be copied and built
        paths :
            src : a list of the paths containing source code. The files will be copied to the build directory (`build.directory`). If a file with the same name is present in several directories in the `build.paths.src` list, then the latest one in the list will overwrite the previous ones.
            remove__ : fortran files named following the pattern <name>__<suffix>.F90 will be renamed in <name>.F90 if this option is set to True (it should probably be!)
        macros : each key should be the name of a *.inc file, and the key values should be the macros that this file defined (e.g. mdf.inc : with_netcdf4 with_go`)
        makefile : name of the makefile to be used
        build_cmd : command to use for building the model (can be within a container!)
        makedepf90 : path to the makedepf90 executable (can be within a container!)

    """
    # Make build directory (and parents, if needed)
    Path(dconf.build.directory).mkdir(exist_ok=True, parents=True)

    if clean :
        [_.unlink() for _ in Path(dconf.build.directory).glob('*.o')]
        [_.unlink() for _ in Path(dconf.build.directory).glob('*.mod')]

    # Copy/Generate source files
    files = copy_files(dconf)
    mfiles = gen_macros(dconf)

    # Remove any source file that would no longer be needed
    cleanup(Path(dconf.build.directory), list(files) + mfiles)

    # Create the makefile
    makefile = gen_makefile(dconf)#.build, dconf.machine.get('host', None))

    # Build TM5
    return make_tm5(dconf) # makefile, dconf.build, dconf.machine.get('host', None))


@debug.trace_call
def copy_files(dconf : DictConfig) -> Dict[str, Path]:
    """
    Compile TM5
    - copy the source files in the build directory
    - create the makefile (with makedepf90)
    - build the executable
    - copy it to the run directory
    :return: path to the compiled executable
    """

    # Ensure that the build directory exists
    build_dir = Path(dconf.build.directory)
    logger.info(f'TM5 will be built in {build_dir}')
    build_dir.mkdir(exist_ok=True)

    # Make a list of all source files to be used:
    files = dict()
    for proj in dconf.build.paths.src:
        for file in Path(proj).iterdir():
            if file.is_file():
                files[file.name] = file

    # Copy the files to the build directory, if they differ from the files that are already in it:
    for filename, source in list(files.items()):
        dest = build_dir / filename

        # Remove the part after '__'? (default True)
        if dconf.build.paths.get("remove__", True):
            if '__' in filename :
                dest = build_dir / (filename.split('__')[0] + '.' + filename.split('.')[-1])
                files[dest.name] = filename

        if not dest.exists():
            logger.info(f"Copy {source} to {dest}")
            shutil.copy(source, dest)
        elif not filecmp.cmp(source, dest):
            logger.info(f"Copy {source} to {dest}")
            shutil.copy(source, dest)

    return files


@debug.trace_call
def gen_macros(dconf: DictConfig) -> List[str]:
    """
    Generate files (*.inc files or sourcecode) that cannot be simply copied from source directories
    """

    # include files:
    fnames = []
    for fname, keys in dconf.build.macros.items():

        # Write the info to a temporary file
        tmpfile = ...
        fid, path = tempfile.mkstemp(text=True)
        with open(path, 'w') as fid:
            for key in keys :
                fid.write(f'#define {key}\n')

        # Copy the file to its final path only if it differs from the existing one:
        dest = Path(dconf.build.directory) / fname
        if not dest.exists():
            shutil.move(fid.name, dest)
        elif not filecmp.cmp(dest, fid.name):
            shutil.move(fid.name, dest)

        # Add the dest file to the list of files regardless
        fnames.append(dest.name)
    return fnames


@debug.trace_call
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


@debug.trace_call
def gen_makefile(dconf: DictConfig) -> Path:
    """
    Create the Makefile_deps file and link to the correct makefile depending on make settings
    """

 #   # Choose the correct Makefile
    buildpath = Path(dconf.build.directory)
    makefile_source = buildpath / dconf.build.makefile
    makefile_dest = buildpath / 'Makefile'
    shutil.copy(makefile_source, makefile_dest)

    # Append it with makedepf90
    # Temporarily move to the build directory for that
    curdir = os.getcwd()
    os.chdir(buildpath)

    command = dconf.build.makedepf90.split() + ['-o', 'tm5.x'] + list(Path().glob('*.[Ff]*'))
    with open(makefile_dest.name, 'a') as fid:  
        subprocess.run(command, stdout=fid)

    os.chdir(curdir)

    return makefile_dest


@debug.trace_call
def make_tm5(dconf: DictConfig) -> Path: #makefile : Path, settings : Dict, host : DictConfig = None) -> Path:
    # Remove the existing executable, to ensure that the script stops if the complilation fails
    executable = Path(dconf.build.directory) / 'tm5.x'
    executable.unlink(missing_ok=True)
    command = dconf.build.build_cmd.split()
    runcmd(command)
    return executable


if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('--build', action='store_true', default=False, help='compile tm5')
    p.add_argument('config', nargs=1)
    args = p.parse_args()

    conf = OmegaConf.load(args.config)
    if args.build :
        build_tm5(conf)

