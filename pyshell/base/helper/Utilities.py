import os, os.path, glob, fnmatch, itertools, colorsys
from datetime import datetime, timedelta
from copy import deepcopy
from numpy import linalg, unique, ma
from numpy import shape, zeros_like, diag, dot
from UserDict import DictMixin
from scipy.interpolate import UnivariateSpline

class ANSI_Colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

    def warn(self, msg, bold=True):
        if bold:
            return self.BOLD + self.WARNING + msg + self.ENDC
        else:
            return self.WARNING + msg + self.ENDC

    def fail(self, msg, bold=True):
        if bold:
            return self.BOLD + self.FAIL + msg + self.ENDC
        else:
            return self.FAIL + msg + self.ENDC

    def emph(self, msg, bold=True):
        if bold:
            return self.BOLD + self.OKBLUE + msg + self.ENDC
        else:
            return self.OKBLUE + msg + self.ENDC

    def ok(self, msg, bold=True):
        if bold:
            return self.BOLD + self.OKGREEN + msg + self.ENDC
        else:
            return self.OKGREEN + msg + self.ENDC

    def uline(self, msg):
        return self.UNDERLINE + msg + self.ENDC

class NearestCorr:
    """
    Computes the nearest correlation matrix to a real symmetric matrix A,
    following the algorithm of Higham et al [2002]
    """
    def __init__(self):
        pass

    def projU(self,U):
        n = shape(U)[0]
        for i in range(n):
            U[i,i] = 1.0
        return U

    def projS(self,A):
        evals, evecs = linalg.eig(A)
        n = shape(evals)[0]
        for i in range(n):
            evals[i] = max(evals[i], 0.0)
        evals = diag(evals)
        M = dot(dot(evecs,evals), evecs.T)
        return M

    def __call__(self, A):
        n = shape(A)[0]
        deltaS = zeros_like(A)
        Y = A
        X = A
        R = A
        for i in range(40):
            R = Y - deltaS
            X = self.projS(R)
            deltaS = X - R
            Y = self.projU(X)
        return Y

def mask_rowcols(a, axis=None):
    """
    Mask rows and/or columns of a 2D array that contain masked values.

    Mask whole rows and/or columns of a 2D array that contain
    masked values.  The masking behavior is selected using the
    `axis` parameter.

      - If `axis` is None, rows *and* columns are masked.
      - If `axis` is 0, only rows are masked.
      - If `axis` is 1 or -1, only columns are masked.

    Parameters
    ----------
    a : array_like, MaskedArray
        The array to mask.  If not a MaskedArray instance (or if no array
        elements are masked).  The result is a MaskedArray with `mask` set
        to `nomask` (False). Must be a 2D array.
    axis : int, optional
        Axis along which to perform the operation. If None, applies to a
        flattened version of the array.

    Returns
    -------
    a : MaskedArray
        A modified version of the input array, masked depending on the value
        of the `axis` parameter.

    Raises
    ------
    NotImplementedError
        If input array `a` is not 2D.

    Notes
    -----
    The input array's mask is modified by this function.

    Examples
    --------
    >>> import numpy.ma as ma
    >>> a = np.zeros((3, 3), dtype=np.int)
    >>> a[1, 1] = 1
    >>> a
    array([[0, 0, 0],
           [0, 1, 0],
           [0, 0, 0]])
    >>> a = ma.masked_equal(a, 1)
    >>> a
    masked_array(data =
     [[0 0 0]
     [0 -- 0]
     [0 0 0]],
          mask =
     [[False False False]
     [False  True False]
     [False False False]],
          fill_value=999999)
    >>> ma.mask_rowcols(a)
    masked_array(data =
     [[0 -- 0]
     [-- -- --]
     [0 -- 0]],
          mask =
     [[False  True False]
     [ True  True  True]
     [False  True False]],
          fill_value=999999)

    """
    a = ma.asarray(a)
    if a.ndim != 2:
        raise NotImplementedError, "compress2d works for 2D arrays only."
    m = ma.getmask(a)
    # Nothing is masked: return a
    if m is ma.nomask or not m.any():
        return a
    maskedval = m.nonzero()
    a._mask = a._mask.copy()
    if not axis:
        a[unique(maskedval[0])] = ma.masked
    if axis in [None, 1, -1]:
        a[:,unique(maskedval[1])] = ma.masked
    return a

def checkDir(fileName, is_dir=False):
    """
    Checks whether the directory corresponding to a file exists.
    If it doesn't, creates it.
    """
    if not is_dir:
        dir_name = os.path.dirname(fileName)
        if not os.path.isdir(dir_name):
            os.makedirs(dir_name)
    else:
        if not os.path.isdir(fileName):
            os.makedirs(fileName)

def readlsmask():
    # read in land/sea mask.
    basemap_datadir = '/home/sbasu/pythonpackages/lib/python2.5/site-packages'
    lsmaskf = open(os.path.join(basemap_datadir,'5minmask.bin'),'rb')
    nlons = 4320; nlats = nlons/2
    delta = 360./float(nlons)
    lsmask = reshape(fromstring(lsmaskf.read(), uint8),(nlats,nlons))
    lsmask_lons = arange(-180,180.,delta)
    lsmask_lats = arange(-90.,90+0.5*delta,delta)
    # add cyclic point in longitude
    lsmask, lsmask_lons = addcyclic(lsmask, lsmask_lons)
    nlons = nlons + 1; nlats = nlats + 1
    # add North Pole point (assumed water)
    tmparr = zeros((nlats,nlons),lsmask.dtype)
    tmparr[0:nlats-1,0:nlons] = lsmask
    lsmask = tmparr
    lsmaskf.close()
    return lsmask_lons, lsmask_lats, lsmask

def loess(x, y, newx, weights=None, alpha=0.25, lmbd=2):
    """
    curve fit using local regression
    ysmooth = loess(x,y,newx,alpha,lambda,robustFlag)
    x,y : data points
    newx,ysmooth : fitted points
    alpha : smoothing, typically 0.25 to 1.0
    lmbd : polynomial order 1 or 2
    """
    # Define weights to be unity if not supplied
    if weights == None:
        weights = ones_like(x)
    # Scale x to [0,1] prevent ill conditioned fitting
    x1 = x[0]
    xr = x.ptp()
    x = (x-x1)/xr
    newx = (newx-x1)/xr
    g = empty_like(newx)
    lmbd = int(round(lmbd)) # force polynm order to be integer
    n = len(x) # number of data points
    q = min(max(floor(alpha*n),lmbd+3),n) # used for weight function width > 3 or so
    # perform a fit for each desired x point
    for ii, x_new in enumerate(newx):
        deltax = abs(x_new - x) # distances from this new point to data
        deltaxsort = sort(deltax) # sorted small to large
        qthdeltax = deltaxsort[q-1] # width of weight function
        arg = minimum(deltax/(qthdeltax*max(alpha,1)),1)
        tricube = (1-abs(arg)**3)**3 # weight function for x distance
        index = where(tricube>0)[0] # select points with nonzero weights
        if len(index) > lmbd:
            p = least2(x[index], y[index], lmbd, weights[index]*tricube[index]) # weighted fit parameters
            newg = polyval(p, x_new) # evaluate fit at this new point
        else:
            newg = average(y[index])+6 # keep same
            raise ValueError('Not enough points')
        g[ii] = newg
    return g

def least2(x,y,n,w=None):
    """"
    p = least2(x,y,n,w) finds the coefficients of a polynomial
    p(x) of degree n that fits the data, p(x(i)) ~= y(i),
    in a weighted least-squares sense with weights w(i).
    """
    if w == None:
        w = ones_like(x)
    # remove data for w=0 to reduce computations and storage
    nzindex = where(w != 0.0)[0]
    x = x[nzindex]
    y = y[nzindex]
    w = w[nzindex]
    nw = len(w)
    pow_array = column_stack([arange(i,i+n) for i in arange(n)])
    LHS = sum([u*(v**pow_array) for u,v in zip(w,x)], axis=0)
    vd_array = vstack([x**i for i in arange(n)])
    RHS = dot(vd_array, y*w)
    p = linalg.solve(LHS, RHS, sym_pos=True, overwrite_a=True, overwrite_b=True)
    return p[::-1]

def savgolay(x, y, n, F, d=0):
    """
    Savitzky-Golay filtering
    Data: (x,y)
    interpolating order n (Must be < n-1)
    Frame length: F (Must be < length(x))
    d: differentiation order >0
    """
    d = int(round(d)) # force differentiation order to be integer
    N = len(x)
    xr = x.ptp()
    x = (x-average(x))/xr # normalise x
    xc = xr**d # correction factor for scaling x
    if F%2 != 0:
        F = F + 1 # F must be odd
    if n > F-1:
        n = F/2 # Interpolating order < window width
    F2 = (F-1)/2 # should be integer
    yf = y # dummy start
    for i in range(F2+1, N-F2+1):
        xloc = x[i-F2:i+F2]
        yloc = y[i-F2:i+F2]
        p = polyfit(xloc,yloc,n)
        if d>0:
            p = polyder(p, d)
        yf[i] = polyval(p, x[i])/xc # center point
    return yf

def gaussian_smooth(x, y, half_window, find_peak=False, peak_type='max'):
    xv = linspace(-half_window, half_window, 2*half_window+1)
    sigma = half_window/3.0
    v = (1./sigma/sqrt(2.*pi)) * exp(-(xv*xv)/2./sigma/sigma)
    smoothed = convolve(y, v, mode='valid')
    if not find_peak:
        return x[half_window:-half_window], smoothed
    else:
        spl = UnivariateSpline(x[half_window:-half_window], smoothed, s=0.0)
        peaks = UnivariateSpline(x[half_window:-half_window], spl(x[half_window:-half_window], 1), s=0.0).roots()
        if len(peaks) == 0:
            return x[half_window:-half_window], smoothed, array([]), array([])
        else:
            spl2 = UnivariateSpline(x[half_window:-half_window], spl(x[half_window:-half_window], 2), s=0.0)
            if peak_type == 'max':
                peaks = filter(lambda q: spl2(q) < 0.0, peaks)
            elif peak_type == 'min':
                peaks = filter(lambda q: spl2(q) > 0.0, peaks)
            return x[half_window:-half_window], smoothed, array(peaks), spl(peaks)

def dirGlob(dir, pattern):
    """ File names matching pattern in directory dir."""
    return glob.glob(os.path.join(dir,pattern))

def OrGlob(dir, *patterns):
    """ File names matching any one of a set of patterns in a directory. """
    return [os.path.join(dir,matchingfile) for matchingfile in filter(lambda x: any([fnmatch.fnmatch(x, p) for p in patterns]), os.listdir(dir))]

def AndGlob(dir, *patterns):
    return [os.path.join(dir,matchingfile) for matchingfile in filter(lambda x: all([fnmatch.fnmatch(x, p) for p in patterns]), os.listdir(dir))]

def dirWalk(args=[], topdir=None, func=dirGlob, nest=False, verbose=False):
    allResults = list()
    # current dir
    if verbose:
            print "*** %s" %topdir
    if topdir is None: topdir = os.getcwd()
    results = func(topdir, *args)
    if verbose:
            print "    %s" % results
    allResults.extend(results)
    # possible sub dirs
    names = [os.path.join(topdir, dir) for dir in os.listdir(topdir)]
    dirs = [n for n in names if os.path.isdir(n)]
    if verbose:
            print "--> %s" % [os.path.basename(d) for d in dirs]
    if len(dirs) > 0:
            for dir in dirs:
                    results = dirWalk(args, dir, func, nest, verbose)
                    if nest:
                            allResults.append(results)
                    else:
                            allResults.extend(results)
    # final allResults
    return allResults

def recurseDirGlob(pattern="*.*", topdir=None, nest=False, verbose=False):
    """
    recurseDirGlob("/home/spir/prog/d0", "*.txt", verbose=True)
    recurseDirGlob("/home/spir/prog/d0", "*.txt")
    recurseDirGlob("/home/spir/prog/d0", "*.txt", nest=True)
    """
    allFilenames = list()
    # current dir
    if verbose:
        print "*** %s" %topdir
    if topdir is None: topdir = os.getcwd()
    filenames = dirGlob(topdir, pattern)
    if verbose:
        for filename in [os.path.basename(d) for d in filenames]:
                print "   %s" %filename
    allFilenames.extend(filenames)
    # possible sub dirs
    names = [os.path.join(topdir, dir) for dir in os.listdir(topdir)]
    dirs = [n for n in names if os.path.isdir(n)]
    if verbose:
        print "--> %s" % [os.path.basename(d) for d in dirs]
    if len(dirs) > 0:
        for dir in dirs:
            filenames = recurseDirGlob(pattern, dir, nest, verbose)
            if nest:
                allFilenames.append(filenames)
            else:
                allFilenames.extend(filenames)
    # final result
    return allFilenames

def IdentityTest(Vals):
    if not Vals:
        return True
    i = iter(Vals)
    first = i.next()
    for item in i:
        if first != item:
            return False
    return True

def Smooth1D(x, window_len, window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t = linspace(-2,2,0.1)
    x = sin(t)+random.randn(len(t))*0.1
    y = smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    """
    if x.ndim != 1:
        raise ValueError, "Smooth1D only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len < 3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window has to be one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s=r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w = ones(window_len,'d')
    else:
        w = getattr(np, window)(window_len)
    y = convolve(w/w.sum(), s, mode='same')
    return y[window_len-1:-window_len+1]

def Smooth2D(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    from scipy import signal
    def gauss_kern(size, sizey=None):
        """ Returns a normalized 2D gauss kernel array for convolutions """
        size = int(size)
        if not sizey:
            sizey = size
        else:
            sizey = int(sizey)
        x, y = mgrid[-size:size+1, -sizey:sizey+1]
        g = exp(-(x**2/float(size) + y**2/float(sizey)))
        return g / g.sum()
    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve(im, g, mode='valid')
    return(improc)

def mdot(*args):
   """Multiply all the arguments using matrix product rules.
   The output is equivalent to multiplying the arguments one by one
   from left to right using dot().
   Precedence can be controlled by creating tuples of arguments,
   for instance mdot(a,((b,c),d)) multiplies a (a*((b*c)*d)).
   Note that this means the output of dot(a,b) and mdot(a,b) will differ if
   a or b is a pure tuple of numbers.
   """
   if len(args)==1:
       return args[0]
   elif len(args)==2:
       return _mdot_r(args[0],args[1])
   else:
       return _mdot_r(args[:-1],args[-1])

def _mdot_r(a,b):
   """Recursive helper for mdot"""
   if type(a) == tuple:
       if len(a)>1:
           a = mdot(*a)
       else:
           a = a[0]
   if type(b) == tuple:
       if len(b)>1:
           b = mdot(*b)
       else:
           b = b[0]
   return dot(a,b)

def Nloop(*lol):
    """
    Converts a list of lists (more generally, a sequence of iterators) to a single iterator that yields a tuple every time an
    element is accessed. For example, if g = Nloop(range(5), ['some', 'holy', 'bugger']), then

    for elem in g: print elem,

    will print

    (0, 'some') (0, 'holy') (0, 'bugger') (1, 'some') (1, 'holy') (1, 'bugger') (2, 'some') (2, 'holy')
    (2, 'bugger') (3, 'some') (3, 'holy') (3, 'bugger') (4, 'some') (4, 'holy') (4, 'bugger')
    """
    l=len(lol)
    dims=map(len,lol)
    totl=reduce(lambda x,y:x*y,dims)
    i=0
    idx=[0]*l
    while totl>i:
        a=i
        i+=1
        for j in range(l-1,-1,-1):
            a,b= divmod(a,dims[j])
            idx[j]=lol[j][b]
        yield tuple(idx)

def Cycle(x):
    """
    Converts a list, or anything that can be converted into an iterator, into a cyclic iterator. For example,

    x = Cycle(['red', 'green', 'blue'])
    for i in range(10):  print x.next(),

    will print

    red green blue red green blue red green blue red
    """
    return itertools.cycle(iter(x))

def Hours(x):
    """
    Takes a datetime() object, extracts the time(), and converts it into hours
    """
    t = x.time()
    return t.hour + t.minute/60. + t.second/3600. + t.microsecond/3600000.

def Year2Datetime(Year):
    """
    Takes a decimal year, for example 2002.44352, and converts it into a datetime object
    """
    Days = 365.242199 * (Year % 1) # Sorry, I'm not worrying about leap years right now
    Year = int(Year)
    return datetime(Year,1,1,0,0,0) - timedelta(days=1) + timedelta(days=Days)

def gmtColormap(fileName,mapName='custom_colors',Segments=512):
    try:
      f = open(fileName)
    except:
      print "file ",fileName, "not found"
      return None
    # Converts a colormap in the GMT file format (.cpt) into something that matplotlib can use
    # Original code by Michale Hearne, available at http://www.mail-archive.com/matplotlib-users@lists.sourceforge.net/msg09547.html
    lines = f.readlines()
    f.close()

    lines = [l.strip() for l in lines if len(l.strip()) > 0]

    x = []
    r = []
    g = []
    b = []
    colorModel = "RGB"
    for l in lines:
        ls = l.split()
        if l[0] == "#":
            if ls[-1] == "HSV":
                colorModel = "HSV"
                continue
            else:
                continue
        if ls[0] == "B" or ls[0] == "F" or ls[0] == "N":
            pass
        else:
            x.append(float(ls[0]))
            r.append(float(ls[1]))
            g.append(float(ls[2]))
            b.append(float(ls[3]))
            xtemp = float(ls[4])
            rtemp = float(ls[5])
            gtemp = float(ls[6])
            btemp = float(ls[7])

    x.append(xtemp)
    r.append(rtemp)
    g.append(gtemp)
    b.append(btemp)

    nTable = len(r)
    x = array( x , float32)
    r = array( r , float32)
    g = array( g , float32)
    b = array( b , float32)
    if colorModel == "HSV":
        for i in range(r.shape[0]):
            rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
            r[i] = rr ; g[i] = gg ; b[i] = bb
    if colorModel == "HSV":
        for i in range(r.shape[0]):
            rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
            r[i] = rr ; g[i] = gg ; b[i] = bb
    if colorModel == "RGB":
        r = r/255.
        g = g/255.
        b = b/255.
    xNorm = (x - x[0])/(x[-1] - x[0])

    red = []
    blue = []
    green = []
    for i in range(len(x)):
        red.append([xNorm[i],r[i],r[i]])
        green.append([xNorm[i],g[i],g[i]])
        blue.append([xNorm[i],b[i],b[i]])
    colorDict = {"red":red, "green":green, "blue":blue}
    return LinearSegmentedColormap(mapName,colorDict,Segments)

def constructColorMaps():
    FileList = glob.glob(os.path.expanduser('~/.ipython/Palettes/*.cpt'))
    mapDict = {}
    for File in FileList:
        mapName = os.path.splitext(os.path.basename(File))[0]
        mapDict[mapName] = gmtColormap(File,mapName)
    return mapDict

def ShowMPLColormaps():
    numCols = 4
    numRows = len(cm.cmap_d)/numCols if len(cm.cmap_d) % numCols == 0 else len(cm.cmap_d)/numCols + 1
    x = linspace(0., 1., 500)
    y = array([0., 1.])
    x, y = meshgrid(x, y)
    figure(figsize=(18, 12))
    for i, (k, v) in enumerate(sorted(cm.cmap_d.items())):
        subplot(numRows, numCols, i+1)
        pcolor(x, cmap=v)
        xticks([])
        yticks([])
        text(x.shape[1]/2.,x.shape[0]/2.,k,ha='center',va='center')
    subplots_adjust(left=0.005,right=0.995,bottom=0.005,top=0.995,wspace=0.02)

def ShowGMTColorMaps():
    colormaps = constructColorMaps()
    numCols = 4
    numRows = len(colormaps)/numCols if len(colormaps) % numCols == 0 else len(colormaps)/numCols + 1
    x = linspace(0., 1., 500)
    y = array([0., 1.])
    x, y = meshgrid(x, y)
    figure(figsize=(18, 12))
    for i, (k, v) in enumerate(sorted(colormaps.items())):
        subplot(numRows, numCols, i+1)
        pcolor(x, cmap=v)
        xticks([])
        yticks([])
        text(x.shape[1]/2.,x.shape[0]/2.,k,ha='center',va='center')
    subplots_adjust(left=0.005,right=0.995,bottom=0.005,top=0.995,wspace=0.05)

class DefaultDict(dict):
    """Dictionary with a default value for unknown keys."""
    def __init__(self, default):
        self.default = default

    def __getitem__(self, key):
        try:
            return self.get(key)
        except:
            return self.setdefault(key, deepcopy(self.default))

    def __copy__(self):
        copy = DefaultDict(self.default)
        copy.update(self)
        return copy

class OrderedDict(dict, DictMixin):

    def __init__(self, *args, **kwds):
        if len(args) > 1:
            raise TypeError('expected at most 1 arguments, got %d' % len(args))
        try:
            self.__end
        except AttributeError:
            self.clear()
        self.update(*args, **kwds)

    def settype(self, value):
        self.default = value

    def clear(self):
        self.__end = end = []
        end += [None, end, end]         # sentinel node for doubly linked list
        self.__map = {}                 # key --> [key, prev, next]
        dict.clear(self)

    def __getitem__(self, key):
        if key in self.keys():
            return self.get(key)
        else:
            self.__setitem__(key, deepcopy(self.default))
            return self.get(key)

    def __setitem__(self, key, value):
        if key not in self:
            end = self.__end
            curr = end[1]
            curr[2] = end[1] = self.__map[key] = [key, curr, end]
        dict.__setitem__(self, key, value)

    def __delitem__(self, key):
        dict.__delitem__(self, key)
        key, prev, next = self.__map.pop(key)
        prev[2] = next
        next[1] = prev

    def __iter__(self):
        end = self.__end
        curr = end[2]
        while curr is not end:
            yield curr[0]
            curr = curr[2]

    def __reversed__(self):
        end = self.__end
        curr = end[1]
        while curr is not end:
            yield curr[0]
            curr = curr[1]

    def popitem(self, last=True):
        if not self:
            raise KeyError('dictionary is empty')
        if last:
            key = reversed(self).next()
        else:
            key = iter(self).next()
        value = self.pop(key)
        return key, value

    def __reduce__(self):
        items = [[k, self[k]] for k in self]
        tmp = self.__map, self.__end
        del self.__map, self.__end
        inst_dict = vars(self).copy()
        self.__map, self.__end = tmp
        if inst_dict:
            return (self.__class__, (items,), inst_dict)
        return self.__class__, (items,)

    def keys(self):
        return list(self)

    setdefault = DictMixin.setdefault
    update = DictMixin.update
    pop = DictMixin.pop
    values = DictMixin.values
    items = DictMixin.items
    iterkeys = DictMixin.iterkeys
    itervalues = DictMixin.itervalues
    iteritems = DictMixin.iteritems

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, self.items())

    def copy(self):
        return self.__class__(self)

    @classmethod
    def fromkeys(cls, iterable, value=None):
        d = cls()
        for key in iterable:
            d[key] = value
        return d

    def __eq__(self, other):
        if isinstance(other, OrderedDict):
            return len(self)==len(other) and self.items() == other.items()
        return dict.__eq__(self, other)

    def __ne__(self, other):
        return not self == other

class PickPoints:

    def __init__(self):
        self.CoordinateList = []
        self.cid = connect('button_press_event', self.PickCoordinate)

    def PickCoordinate(self, event):
        x, y = event.x, event.y
        if event.inaxes and event.button == 3:
            self.CoordinateList.append([event.xdata, event.ydata])

    def __call__(self):
        disconnect(self.cid)
        return array(self.CoordinateList)
