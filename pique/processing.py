#!/usr/bin/env python
"""
Data processing functions for Pique.
"""
import numpy
import scipy.signal
import itertools
import operator

# numpy burbles about silly things.
numpy.seterr( all='ignore' )

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

    data = scipy.signal.wiener( data, mysize=alpha )
    data = scipy.signal.convolve( data, window, mode='same' )
    
    return data

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

