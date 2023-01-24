#!/usr/bin/env python

import os
import setuptools


setuptools.setup(
    name="tmflex",
    version="2023.01.16",
    author="Guillaume Monteil",
    author_email="guillaume.monteil@nateko.lu.se",
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
    install_requires=['numpy', 'ipython', 'pandas', 'xarray', 'h5py', 'netcdf4', 'loguru', 'omegaconf', 'mkdocs'],
    scripts=['bin/tm5'],
)
