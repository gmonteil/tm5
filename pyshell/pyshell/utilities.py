#!/usr/bin/env python2.7

import os
from netCDF4 import Dataset


def checkDir(filename, is_dir=False):
    """
    Checks whether the directory corresponding to a file exists.
    If it doesn't, creates it.
    """
    if not is_dir:
        dir_name = os.path.dirname(filename)
        if not os.path.isdir(dir_name):
            os.makedirs(dir_name)
    else:
        if not os.path.isdir(filename):
            os.makedirs(filename)


class my_Dataset(Dataset):
    def __init__(self, file_name, *args, **kwargs):

        # We need to check if the mode is 'r' or 'a'. Basically, if it's not 'w', check for file existence.
        check_exist = True
        if 'mode' in kwargs:
            check_exist = not kwargs['mode'].startswith('w')
        if len(args) > 0:
            check_exist = not args[0].startswith('w')
        if check_exist:
            if not os.path.exists(file_name):
                raise RuntimeError('File %s does not exist'%file_name)
        Dataset.__init__(self, file_name, *args, **kwargs)