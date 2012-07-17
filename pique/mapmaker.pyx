import numpy
cimport numpy
cimport cython
from libc.stdlib cimport *

cdef extern from 'math.h':
    double ceil(double x)
    double floor (double x)

DTYPE = numpy.int # data type for arrays
ctypedef numpy.int_t DTYPE_t

@cython.boundscheck(False)
def hist( numpy.ndarray[DTYPE_t, ndim=1] data,  \
        int lowest,                             \
        int highest,                            \
        int bins,                               \
        int window,                             \
        int stride ) :
    cdef int overlap  = int(ceil( window / float(stride) ))
    cdef int w = len(data)/stride + overlap
    cdef numpy.ndarray r = numpy.zeros( [ bins, w ], dtype=DTYPE )
    cdef int interval = int(ceil( ( highest - lowest ) / float(bins) ))
    cdef int x
    cdef int i
    cdef int j
    cdef int n
    cdef int m
    for x in range(len(data)) :
        i = data[x]
        if i >= lowest and i <= highest :
            j = (i-lowest)/interval
            n = x/stride
            r[ j, n:(n+overlap) ] += 1
    return r
