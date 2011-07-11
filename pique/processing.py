#!/usr/bin/env python
"""
Data processing functions for Pique.
"""
import numpy
import scipy
import itertools
import operator

def findregions( data, N ) :
    """
    Return all the regions of an array that exceed N.
    """
    pos = filter(lambda(a) : a[1] == abs(a[1]), zip(range(len(data)),data-N) )
    if len(pos) == 0 :
        return []
    pos = numpy.array( zip(*pos)[0] )
    regions = []
    for k,g in itertools.groupby( enumerate(pos), lambda(i,x):i-x ) :
        l = map( operator.itemgetter(1),g )
        regions.append( {'start':l[0],'stop':l[-1] } )
    return regions

def sortbyintegral( data, regions ) :
    """
    Sort regions of an array according to their definite integrals.
    """
    return sorted( regions, key=lambda r: sum(data[r['start']:r['stop']]) )
    
def filterset( data, alpha, l_thresh ) :
    """
    Core data filtering function.
    """
    
    window = scipy.signal.blackman( l_thresh )
    window = window / sum(window)

    dataf = scipy.signal.detrend( data )
    dataf = scipy.signal.wiener( dataf, mysize=alpha )
    datar = scipy.signal.detrend( data[::-1] )
    datar = scipy.signal.wiener( datar, mysize=alpha )

    dataf = scipy.signal.convolve( dataf, window, mode='full' )
    datar = scipy.signal.convolve( datar, window, mode='full' )

    return ( dataf + datar[::-1] ) / 2.0

def overlaps( forward_regions, reverse_regions ) :
    """
    Find enriched regions that meet the overlap criterion.
    """
    envelopes = []
    for l in forward_regions :
        for m in reverse_regions :
            if m['start'] < l['stop'] < m['stop'] :
                if l['start'] < m['start'] < l['stop'] :
                    e = {   'forward':l,        \
                            'reverse':m,        \
                            'start':l['start'], \
                            'stop':m['stop'],   \
                            'annotations':{}    }
                    envelopes.append( e )
    return envelopes

