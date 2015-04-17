from libc.stdlib cimport malloc, free
import collections
import numpy as np
import cython
cimport numpy as np

np.import_array()

from stepstr import*

# import c fuctions
cdef extern from "cdtw.c":

    cdef int extra_size(int dp_type)

    cdef struct t_path_element:
        int i
        int j

    cdef struct t_weights:
        double a
        double b
        double c

    cdef struct t_dtw_settings:
        int compute_path
        int dist_type
        int dp_type
        int window_type
        double window_param
        int norm
        int offset
        t_weights weights

    cdef double cdtw(double * ref,
                     double * query,
                     int len_ref,
                     int len_query,
                     double * cost_matrix,
                     t_path_element * path,
                     int * true_path_len,
                     t_dtw_settings dtw_settings)

# macros
cdef enum:
    _EUCLID = 11
    _EUCLID_SQUARED = 12
    _MANHATTAN = 13
    _DPW = 20
    _DP1 = 21
    _DP2 = 22
    _DP3 = 23
    _SCP0SYM = 24
    _SCP0ASYM = 25
    _SCP1DIV2SYM = 26
    _SCP1DIV2ASYM = 27
    _SCP1SYM = 28
    _SCP1ASYM = 29
    _SCP2SYM = 210
    _SCP2ASYM = 211

    _SCBAND = 31
    _PALIVAL = 32
    _ITAKURA = 33
    _PALIVAL_MOD = 35


cdef path_wrapper(t_path_element * cpath, int cpath_len):
    """Path wrapper

    Converts path C array to python list of tuples.

    Args:
        t_path_element *cpath: pointer to t_path_elements array
        int cpath_len: length og the cpath array

    Returns:
        list: list of tuples, each tuple contains row and columns indices    

    """
    path = []
    cdef int m
    for m in range(cpath_len):
        path.append((cpath[cpath_len - m - 1].i, cpath[cpath_len - m - 1].j))
    return path


class Dist:

    """Distance type class

    Contains dintance types for dynamic time wapring algorithm. There are three 
    available distance functions at the moment: 'manhattan', 'euclid' and 
    'euclid_squared'.
    """

    def __init__(self, dist='euclid'):
        """__init__ method

        Args:

            dist(str, optional): distance type, default is 'manhattan'

        """
        self._cur_type = dist
        self._types = {'euclid': _EUCLID,
                       'euclid_squared': _EUCLID_SQUARED,
                       'manhattan': _MANHATTAN}

    def set_type(self, itype):
        """set_type method

        Set distance type, available distances are: 
        'manhattan', 'euclid' and 'euclid_squared'

        Args:

            itype(str): distance type to set

        """
        if itype in self._types.keys():
            self._cur_type = itype
        else:
            print "Unknown distance function type, \
                   available distance functions\n" + str(self._types)

    def get_cur_type(self):
        return self._types[self._cur_type]

    def __str__(self):
        return str(self._cur_type)


class Step:

    """Step class

    Class containts different step patterns for dynamic time warping algorithm.
    There are folowing step patterns available at the moment:
    Usual step patterns:
        'dp1' 
        'dp2' 
        'dp3'
    Sakoe-Chiba classification:
        'p0sym':
        'p0asym':
        'p05sym':
        'p05asym':
        'p1sym':
        'p1asym':
        'p2sym':
        'p2asym':

    You can see step pattern definition using print_step method
    """

    def __init__(self, step='dp2', weights = [1,1,1]):
        self._cur_type = step
        self._types = { 'dpw': _DPW, 'dp1': _DP1, 'dp2': _DP2, 'dp3': _DP3,
                       'p0sym': _SCP0SYM, 'p0asym': _SCP0ASYM,
                       'p05sym': _SCP1DIV2SYM, 'p05asym': _SCP1DIV2ASYM,
                       'p1sym':  _SCP1SYM, 'p1asym': _SCP1ASYM,
                       'p2sym': _SCP2SYM, 'p2asym': _SCP2ASYM}
        self._weights = weights

    def set_type(self, itype):
        if itype in self._types.keys():
            self._cur_type = itype
        else:
            print "Unknown DP type, available slope constraints: \n"  \
                + str(self._types)

    def set_weights(self,w):
        self._weights = w

    def get_cur_type(self):
        return self._types[self._cur_type]

    def get_weights(self):
        return self._weights

    def step_str(self, itype):
        return stepstr[itype]

    def __str__(self):
        return str(self._cur_type)


class Window:

    """Global constraint class.
    Available constraints: scband, itakura, palival, itakura_mod
    """

    def __init__(self, window='nowindow', param=0.0):
        self._cur_type = window
        self._types = {'scband': _SCBAND, 'palival':     _PALIVAL,
                       'itakura': _ITAKURA, 'palival_mod': _PALIVAL_MOD, 'nowindow': 0}
        self._dummy = {_SCBAND:'scband', _PALIVAL: 'palival',\
                       _ITAKURA:'itakura', _PALIVAL_MOD:'palival_mod', 0: 'nowindow'}
        self._param = param

    def set_type(self, itype):
        if itype in self._types.keys():
            self._cur_type = itype
        else:
            print "Unknown Global Constraint type, \
            available slope constraints: \n" + str(self.__types.keys())

    def set_param(self, param):
        self._param = float(param)

    def get_cur_type(self):
        return self._types[self._cur_type]

    def get_param(self):
        return self._param

    def get_options(self):
        return self._types

    def __str__(self):
        return str(self._cur_type + ', parameter: ' + str(self._param))


class Settings:

    """
    class with dtw settings
    """

    def __init__(self,
                 dist='manhattan',
                 step='dp2',
                 window='nowindow',
                 param=0.0,
                 norm=False,
                 compute_path=False):

        self.dist = Dist(dist)
        self.step = Step(step)
        self.window = Window(window, param)
        self.compute_path = compute_path
        self.norm = norm

    def __str__(self):
        return str('distance function: ' + str(self.dist) + '\n'
                   'local constraint: ' + str(self.step) + '\n'
                   'window: ' + str(self.window) + '\n'
                   'normalization: ' + str(self.norm) + '\n')


class cydtw:

    """
    Main entry, dtw algorithm
    settings is optional parameter
    default settings are dp2 without global constraint
    """

    def __init__(self, ref, query, settings=Settings()):
        self._dist = None
        self._cost = [[]]
        self._dir = [[()]]
        self._path = [()]
        self._dtw(ref, query, settings)

    def _dtw(self, ref, query, settings):
        # sequence control
        if len(ref) == 0 or len(query) == 0:
            return
        # map python settings to C structure dtw_settings functions
        cdef t_dtw_settings c_dtw_settings
        c_dtw_settings.compute_path = <int > settings.compute_path
        c_dtw_settings.dist_type = <int > settings.dist.get_cur_type()
        c_dtw_settings.dp_type = <int > settings.step.get_cur_type()
        c_dtw_settings.window_type = <int > settings.window.get_cur_type()
        c_dtw_settings.window_param = <double > settings.window.get_param()
        c_dtw_settings.norm = <int > settings.norm
        c_dtw_settings.offset = extra_size(c_dtw_settings.dp_type)
        c_dtw_settings.weights.a = <double>settings.step.get_weights()[0]
        c_dtw_settings.weights.b = <double>settings.step.get_weights()[1]
        c_dtw_settings.weights.c = <double>settings.step.get_weights()[2]

        # allocate path
        cdef t_path_element * cpath = <t_path_element * >malloc(sizeof(
                               t_path_element)*(< int > ( len(ref) + len(query))))

        # init true path length
        cdef int cpath_len = 0

        # expand ref and query
        if isinstance(ref, np.ndarray) and isinstance(query, np.ndarray):
            ref = np.hstack((np.zeros((c_dtw_settings.offset)), ref))
            query = np.hstack((np.zeros((c_dtw_settings.offset)), query))
        else:
            ref = [0 for _ in range(c_dtw_settings.offset)] + ref
            query = [0 for _ in range(c_dtw_settings.offset)] + query

        # init numpy arrays
        cdef np.ndarray[np.float_t, ndim = 1] cref
        cdef np.ndarray[np.float_t, ndim = 1] cquery
        cdef np.ndarray[np.float_t, ndim = 2] cost

        # contiguous c array in memory
        cref = np.ascontiguousarray(ref, dtype=np.float)
        cquery = np.ascontiguousarray(query, dtype=np.float)

        # init cost matrix
        cost = np.zeros((cref.shape[0], cquery.shape[0]),
                        dtype=np.float)
        # call cdtw function (in cdtw.c)
        self._dist =  cdtw( < double*>cref.data,
                           < double * >cquery.data,
                           < int > len(ref) - c_dtw_settings.offset,
                           < int > len(query) - c_dtw_settings.offset,
                           < double*>&cost[0, 0],
                           cpath,
                           & cpath_len,
                           c_dtw_settings);

        self._cost = cost[c_dtw_settings.offset:, c_dtw_settings.offset:]

        # convert c path to python path
        if(settings.compute_path):
            self._path = path_wrapper(cpath, cpath_len)

        # cleaning
        free(cpath)

    def get_dist(self):
        return self._dist

    def get_cost(self):
        return self._cost

    def get_path(self):
        return self._path
