#!/usr/bin/env python
"""
Generate coverage tracks from a BAM file.
"""
import pysam
import numpy

def loadBAM( file ) :
    samfile = pysam.Samfile( file, 'rb' )
    tracks = {}
    for replicon in samfile.header['SQ'] :
        name   = replicon['SN']
        length = replicon['LN']
        tracks[name] = {    'length'  : length,                 \
                            'forward' : numpy.zeros(length),    \
                            'reverse' : numpy.zeros(length), }
        for pileupcolumn in samfile.pileup( name, 0, length) :
            fn,rn = 0,0
            for pileupread in pileupcolumn.pileups :
                if pileupread.alignment.is_reverse :
                    rn = rn + 1
                else :
                    fn = fn + 1
            tracks[name]['forward'][pileupcolumn.pos] = fn
            tracks[name]['reverse'][pileupcolumn.pos] = rn
    return tracks
