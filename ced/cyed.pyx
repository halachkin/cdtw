cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free

np.import_array()

from stepstr import *

# import c fuctions
cdef extern from "ced.c":
    cdef int extra_size(int dp_type)

    cdef struct t_path_element:
        int i
        int j

    cdef struct t_weights:
        double a
        double b
        double c

    cdef struct t_extra_args:
        double * sigmas
        double sigma
        double * gap

    cdef struct t_settings:
        int compute_path
        int dist_type
        int qtse_type
        int dp_type
        int window_type
        double window_param
        int norm_type
        int offset
        t_weights weights

    cdef double ced(double *ref,
                    double *query,
                    t_extra_args args,
                    int len_ref,
                    int len_query,
                    int ncols,
                    double *cost_matrix,
                    t_path_element *path,
                    int *true_path_len,
                    t_settings settings)

# macros
cdef enum:
    _EUCLID = 11
    _EUCLID_SQUARED = 12
    _MANHATTAN = 13

    _DPW = 20
    _DP1 = 21
    _DP2 = 22
    _DP2_EDR = 220
    _DP2_ERP = 221
    _DP3 = 23
    _DP3_LCSS = 230
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

    _NO_QUANTISATION = 40
    _EDR = 41
    _LCSS = 42

    _NO_NORMALISATION = 50
    _NORM_BY_MIN_LENGTH = 51
    _NORM_BY_AVG_LENGTH = 52
    _NORM_BY_MAX_LENGTH = 53


cdef path_wrapper(t_path_element *cpath, int cpath_len):
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


class Setting:
    def __init__(self):
        self._types = {}
        self._cur_type = str()

    def set_type(self, itype):
        """set type method"""
        if itype in self._types.keys():
            self._cur_type = itype
        else:
            print "Unknown type, possible types: \n" + str(self._types.keys())

    def get_cur_type(self):
        """get type method"""
        return self._cur_type

    def get_cur_type_code(self):
        return self._types[self._cur_type]

    def get_options(self):
        """get options method"""
        return self._types.keys()

    def __str__(self):
        """__str__ method"""
        return str(self._cur_type)


class Dist(Setting):
    """Distance type class
    Contains dintance types for dynamic time wapring algorithm. There are three 
    available distance functions at the moment: 
    'manhattan', 
    'euclid',
    'euclid_squared'.
    """

    def __init__(self, dist='euclid'):
        """__init__ method
        Args:
            dist(str, optional): distance type, default is 'manhattan'
        """
        Setting.__init__(self)
        self._cur_type = dist
        self._types = {'euclid': _EUCLID,
                       'euclid_squared': _EUCLID_SQUARED,
                       'manhattan': _MANHATTAN}


class Quantisation(Setting):
    def __init__(self, qtse='no_quantisation'):
        """__init__ method
        Args:
            qtse(str, optional): quantisation type, default is 'no_quantisation'
        """
        Setting.__init__(self)
        self._cur_type = qtse
        self._types = {'edr': _EDR, 'lcss': _LCSS,
                       'no_quantisation': _NO_QUANTISATION}


class Normalisation(Setting):
    def __init__(self, norm='none'):
        """__init__ method
        Args:
            qtse(str, optional): quantisation type, default is 'no_quantisation'
        """
        Setting.__init__(self)
        self._cur_type = norm
        self._types = {'none': _NO_NORMALISATION, 'min': _NORM_BY_MIN_LENGTH,
                       'avg': _NORM_BY_AVG_LENGTH, 'max': _NORM_BY_MAX_LENGTH}


class Step(Setting):
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

    def __init__(self, step='dp2', weights = [1, 1, 1]):
        self._cur_type = step
        self._types = {'dpw': _DPW, 'dp1': _DP1, 'dp2': _DP2, 'dp2edr': _DP2_EDR, 'dp2erp': _DP2_ERP,
                       'dp3': _DP3, 'dp3_lcss': _DP3_LCSS,
                       'p0sym': _SCP0SYM, 'p0asym': _SCP0ASYM,
                       'p05sym': _SCP1DIV2SYM, 'p05asym': _SCP1DIV2ASYM,
                       'p1sym': _SCP1SYM, 'p1asym': _SCP1ASYM,
                       'p2sym': _SCP2SYM, 'p2asym': _SCP2ASYM}
        self._weights = weights

    def set_weights(self, w):
        self._weights = w

    def get_weights(self):
        return self._weights

    def step_str(self, itype):
        return stepstr[itype]

    def __str__(self):
        return str(self._cur_type)


class Window(Setting):
    """Global constraint class.
    Available constraints: scband, itakura, palival, itakura_mod
    """

    def __init__(self, window='nowindow', param=0.0):
        self._cur_type = window
        self._types = {'scband': _SCBAND, 'palival': _PALIVAL,
                       'itakura': _ITAKURA, 'palival_mod': _PALIVAL_MOD, 'nowindow': 0}

        self._param = param

    def set_param(self, param):
        self._param = float(param)

    def get_param(self):
        return self._param

    def __str__(self):
        return str(self._cur_type + ', parameter: ' + str(self._param))


class Settings:
    """
    class with ed settings
    """

    def __init__(self,
                 dist='manhattan',
                 qtse='no_quantisation',
                 step='dp2',
                 window='nowindow',
                 param=0.0,
                 norm='no_normalisation',
                 compute_path=False):
        self.dist = Dist(dist)
        self.step = Step(step)
        self.qtse = Quantisation(qtse)
        self.window = Window(window, param)
        self.compute_path = compute_path
        self.norm = Normalisation(norm)

    def __str__(self):
        return 'Distance function: {}\nLocal constraint: {}\nWindow: {}\nNormalisation: {}\n'\
            .format(self.dist, self.step, self.window, self.norm)


class cyed:
    """
    Main entry, ed algorithm
    settings is optional parameter
    default settings are dp2 without global constraint
    """

    def __init__(self, ref, query, args, settings=Settings()):
        self._dist = None
        self._cost = [[]]
        self._dir = [[()]]
        self._path = [()]
        self._ed(ref, query, args, settings)

    @staticmethod
    def flatten_padded(arr, offset):
        arr_shape = np.shape(arr)
        ndims = len(arr_shape)
        if ndims > 2:
            raise Exception('Multi-dimensional input must be a two-dimensional matrix')

        if ndims == 1:
            pad_size = (offset, )
            return np.hstack((np.zeros(pad_size), arr)).ravel(order='C')
        else:
            pad_size = (offset, arr_shape[1])
            return np.vstack((np.zeros(pad_size), arr)).ravel(order='C')


    def _ed(self, ref, query, args, settings):
        if not isinstance(ref, np.ndarray) or not isinstance(query, np.ndarray):
            raise Exception('Sequences must be numpy arrays')

        # sequence control
        if len(ref) == 0 or len(query) == 0:
            raise Exception('Length of both sequences must be > 0')

        ref_shape = np.shape(ref)
        qry_shape = np.shape(query)

        if len(ref_shape) != len(qry_shape):
            raise Exception('Two sequences must have the same number of dimensions')

        ndims = len(ref_shape)

        if ndims > 2:
            raise Exception('Currently only supports one or two dimensional sequences')

        if ndims == 2 and ref_shape[1] != qry_shape[1]:
            raise Exception('Shapes of two sequence can only differ in the first dimension')

        # map python settings to C structure settings functions
        cdef t_settings c_settings
        c_settings.compute_path = <int> settings.compute_path
        c_settings.dist_type = <int> settings.dist.get_cur_type_code()
        c_settings.qtse_type = <int> settings.qtse.get_cur_type_code()
        c_settings.dp_type = <int> settings.step.get_cur_type_code()
        c_settings.window_type = <int> settings.window.get_cur_type_code()
        c_settings.window_param = <double> settings.window.get_param()
        c_settings.norm_type = <int> settings.norm.get_cur_type_code()
        c_settings.offset = extra_size(c_settings.dp_type)
        c_settings.weights.a = <double> settings.step.get_weights()[0]
        c_settings.weights.b = <double> settings.step.get_weights()[1]
        c_settings.weights.c = <double> settings.step.get_weights()[2]

        # allocate path
        cdef t_path_element *cpath = <t_path_element *> malloc(sizeof(
                                                               t_path_element) * (<int> (len(ref) + len(query))))

        # init true path length
        cdef int cpath_len = 0
        len_ref = ref_shape[0]
        len_qry = qry_shape[0]
        offset = c_settings.offset

        if ndims == 1:
            ncols = 1
        else:
            ncols = ref_shape[1]

        # init numpy arrays
        cdef np.ndarray[np.float_t, ndim = 1] cref
        cdef np.ndarray[np.float_t, ndim = 1] cquery
        cdef np.ndarray[np.float_t, ndim = 1] csigmas
        cdef np.ndarray[np.float_t, ndim = 1] cgap
        cdef np.ndarray[np.float_t, ndim = 2] cost

        # expand ref and query and ravel as one dimensional, contiguous c arrays in memory
        cref = self.flatten_padded(ref, offset)
        cquery = self.flatten_padded(query, offset)

        cdef t_extra_args cargs
        if 'sigmas' in args:
            csigmas = self.flatten_padded(args['sigmas'], 0)
            cargs.sigmas = <double*> csigmas.data

        if 'sigma' in args:
            cargs.sigma = <double> args['sigma']

        if 'gap' in args:
            cgap = self.flatten_padded(args['gap'], 0)
            cargs.gap = <double*> cgap.data

        # init cost matrix
        cost = np.zeros((len_ref + offset, len_qry + offset),
                        dtype=np.float)
        # call ced function (in ced.c)
        self._dist = ced(<double*> cref.data,
                         <double *> cquery.data,
                         cargs,
                         <int> len_ref,
                          <int> len_qry,
                          <int> ncols,
                          <double*> &cost[0, 0],
                          cpath,
                          &cpath_len,
                          c_settings)

        self._cost = cost[offset:, offset:]

        # convert c path to python path
        if (settings.compute_path):
            self._path = path_wrapper(cpath, cpath_len)

        # cleaning
        free(cpath)

    def get_dist(self):
        return self._dist

    def get_cost(self):
        return self._cost

    def get_path(self):
        return self._path
