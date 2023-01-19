#!/usr/bin/env python

import h5py
from numpy import array, append, zeros, unique
from datetime import datetime


class observations:
    def __init__(self):
        self.n = 0
        self.tracers = []
        self.data = dict.fromkeys(['latitudes', 'longitudes', 'altitudes', 'mixing_ratio', 'mixing_ratio_err', 'sampling_strategy', 'dates', 'sitecode', 'tracer', 'dates'], array(()))

    def add(self, tracer, pattern, obsclass):
        import glob
        flist = glob.glob(pattern)
        for ff in flist :
            self.append(obsclass(ff))
        if tracer not in self.tracers:
            self.tracers.append(tracer)

    def append(self, obsobj):
        for key in self.data :
            try :
                self.data[key] = append(self.data[key], obsobj.data[key])
            except KeyError :
                if key == 'tracer' :
                    self.data['tracer'] = append(self.data['tracer'], array([obsobj.param]*obsobj.n))
                else :
                    import pdb; pdb.set_trace()
        self.n += obsobj.n

    def get(self, tracer, ti, tf):
        outdict = {}
        # add "station_id" field needed by TM5. Do this before the filtering by date, so that in case of input file splitting, the indexing remains consistent TODO: replace this by a proper station database.
        statid = zeros(self.n)
        for istat, stat in enumerate(unique(self.data['sitecode'])) :
            statid[self.data['sitecode'] == stat] = istat
        self.data['tracer'] = array([x.lower() for x in self.data['tracer']])   # make sure the filter won't be case sensitive

        # select the obs that are in the requested time period
        selection = self.data['tracer'] == tracer.lower()
        selection = selection*(self.data['dates'] >= ti)*(self.data['dates'] <= tf)

        # fill in the dictionary to pass on to TM5:
        for key in self.data :
            try :
                outdict[key] = self.data[key][selection]
            except :
                import pdb; pdb.set_trace()
        outdict['station_id'] = statid[selection]
        outdict['n'] = sum(selection)
        return outdict


class HDF5DB:
    def __init__(self, filename):
        df = h5py.File(filename, 'r')
        self.n = df.attrs['n']
        self.param = 'co2'
        self.data = {
                'latitudes': df['lat'][:],
                'longitudes': df['lon'][:],
                'altitudes': df['alt'][:],
                'mixing_ratio': df['obs'][:],
                'mixing_ratio_err': df['err'][:],
                'dates': array([datetime(*x) for x in df['time']]),
                'sitecode': df['site'][:],
                'tracer': df['tracer'][:],
                'sampling_strategy': df['sampling_strategy'][:]
                }
        if 'obspack_id' in df :
            self.data['obspack_id'] = df['obspack_id'][:]