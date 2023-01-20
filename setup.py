#!/usr/bin/env python

import os
import setuptools


setuptools.setup(
    name="pyshell",
    version="2023.01.16",
    author="Guillaume Monteil",
    author_email="guillaume.monteil@nateko.lu.se",
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2",
        "Operating System :: OS Independent",
    ],
    python_requires='>=2.7',
    install_requires=['h5py', 'scipy', 'python-dateutil', 'netcdf4', 'progressbar'],
    extras_require={'interactive': ['ipython']},
    scripts=['bin/pyshell', 'bin/submit_tm5_step_run', 'bin/submit_tm5', 'bin/setup_tm5', 'bin/congrad.exe'],
)
