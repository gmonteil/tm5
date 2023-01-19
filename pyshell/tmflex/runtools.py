#!/bin/env python

from pyshell.tmflex.rc import RcFile
import os, re

def parse_args():
    parser = ArgumentParser(usage='Usage: %(prog)s <options> action rc-file', formatter_class=lambda prog: HelpFormatter(prog,max_help_position=30,width=120))
    
    parser.add_argument('action', choices=['initiate', 'setup', 'forward', 'adjoint', 'finish', 'restart'], metavar='action')
    args = parser.parse_args()
    return args

class rcdat(RcFile):
    def __init__(self, filename=None, silent=False, marks=('${','}')):
        RcFile.__init__(self, filename, silent, marks)
        self.replace_add = self.setkey
        self.filename=filename

    def substituteTimes(self, StartTime=None, EndTime=None):
        """
        The rc file is allowed to have character combinations such as <Y1> and <m2>, which mean,
        respectively, the starting year and the ending month. This routine goes through all the
        rc keys and makes the necessary substitutions.
        """
        if StartTime == None :
            try :
                StartTime = self.ti
            except AttributeError :
                raise
        if EndTime == None :
            try :
                EndTime = self.tf
            except AttributeError :
                raise
        match_string = []
        for atom in ['Y', 'y', 'm', 'M', 'd', 'H', 'b', 'B', 'j', 'p', 'S']:
            for i in range(1,3):
                rep_string = '<%1s%1i>'%(atom, i)
                match_string.append(rep_string)
        match_string = '|'.join(match_string)
        match_pattern = re.compile(match_string)
        for key in self.keys():
            value = self.get(key)
            if type(value) == str :
                if match_pattern.search(value):
                    value = putDateString(StartTime, EndTime, value)
                    self.replace(key, value)
    def copy(self):
        from copy import deepcopy
        return deepcopy(self)
    def setkey(self, key, val):
        if self.has_key(key):
            self.replace(key, val)
        else :
            self.add(key, val)
    def save(self, path='./rcf.pic'):
        import pickle
        file_name = path
        with open(file_name, 'wb') as fid:
            pickle.dump(self, fid, pickle.HIGHEST_PROTOCOL)

    def setup_meteo_coarsening(self, coarsen_meteo):
        if coarsen_meteo:
            self.setkey('my.meteo.resol', 'glb100x100')
            self.setkey('my.tmm.output', True)
        else :
            self.setkey('my.meteo.resol', 'coarsened')
            self.setkey('my.tmm.output', False)

class archive:
    def __init__(self, host, path, protocol='rsync', options=''):
	self.host = host    # remote host
	self.rpath = path   # remote path
	self.protocol = protocol # protocol: ssh/scp or rsync
	self.paths = []     # paths to archive
	self.names = []     # optional renaming of the paths
	self.options = options

    def add(self, path, name=None):
	self.paths.append(path)
	if name == None : name = os.path.split(path)[1]
	self.names.append(name)

    def save(self):
	for path, name in zip(self.paths, self.names):
	    if self.protocol == 'rsync' :
		os.system('rsync %s -avh %s %s:%s/%s'%(self.options, path, self.host, self.rpath, name))
	    elif self.protocol in ['scp', 'ssh'] :
		os.system('scp %s -r %s %s:%s/%s'%(self.options, path, self.host, self.rpath, name))

    def cleanup(self):
        self.paths = []
        self.names = []

def putDateString(StartTime, EndTime, input_string):
    output_string = input_string
    for atom in ['Y', 'y', 'm', 'M', 'd', 'H', 'b', 'B', 'j', 'p', 'S']:
	for i, t in zip([1,2], [StartTime, EndTime]):
	    rep_string = '<%1s%1i>'%(atom, i)
	    form_string = '%%%1s'%atom
	    output_string = output_string.replace(rep_string, t.strftime(form_string))
    return output_string

