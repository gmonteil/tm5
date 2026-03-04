#!/usr/bin/env python
### collection of general utilities
###

import numpy as np
import xxhash


def md5(fname: str, chunk_size: int=1024 * 1024):
    h = xxhash.xxh32()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()

def utc_to_lst( time_utc : np.ndarray | np.datetime64, longitude : np.ndarray | float ):
    """conversion of UTC time to local solar time (LST),
    longitude coordinate may be provided as
    scalar (equally applied for all time points)
    array  (individual longitude for each time point)
    """
    #-- ensure same length (if both are arrays)
    if hasattr(time_utc, "__len__") and hasattr(longitude, "__len"):
        assert len(time_utc)==len(longitude)
        
    #-- UTC/LST time difference (15 degrees equals 1 hour)
    if hasattr(longitude, "__len__"):
        td = np.asarray(longitude/15. * 3600, dtype='i4') # [s]
        #-- as timedelta64
        td = np.array(td, dtype='timedelta64[s]')
    else:
        td = int(longitude/15. * 3600) # [s]
        td = np.timedelta64(td, 's')
    #--
    time_LST = time_utc + td

    return time_LST
