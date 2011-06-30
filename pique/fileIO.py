#!/usr/bin/env python
"""
Generate coverage tracks from a BAM file.
"""
import pysam
import numpy
import sys

GFFkeys = ['contig','source','feature','start','stop','score','strand','frame','group']


def loadBAM( file ) :
    """
    Read eeach contig in a BAM file and return a dictionary of tracks.
    Each track contains scalar length, and two vectors representing
    the forward and everse coverage across the contig.
    """
    samfile = pysam.Samfile( file, 'rb' )
    tracks = {}

    # loop over the congigs and build forward and reverse coverage
    # tracks for them
    for n,contig in enumerate(samfile.references) :
        
        length = samfile.lengths[n]
        tracks[contig] = {  'length'  : length,                 \
                            'forward' : numpy.zeros(length),    \
                            'reverse' : numpy.zeros(length), }
        
        for pileupcolumn in samfile.pileup( contig, 0, length) :
            
            fn,rn = 0,0
            
            for pileupread in pileupcolumn.pileups :
                
                if pileupread.alignment.is_reverse :
                    rn = rn + 1
                else :
                    fn = fn + 1
                
            tracks[contig]['forward'][pileupcolumn.pos] = fn
            tracks[contig]['reverse'][pileupcolumn.pos] = rn
            
    samfile.close()
    
    return tracks

def loadGFF( file ) :
    """
    Parse a GFF file and return the analysis region and mask features.
    """
    regions = []
    masks   = []
    for line in open( file ) :
       g = dict(zip( GFFkeys, line.strip().split('\t') ))
       if g['feature'] == 'analysis_region' :
            region = {  'contig' : g['contig']      ,   \
                        'start'  : int( g['start'] ),   \
                        'stop'   : int( g['stop']  )    }
            regions.append(region)

       if g['feature'] == 'mask' :
            mask =  {   'contig' : g['contig']      ,   \
                        'start'  : int( g['start'] ),   \
                        'stop'   : int(g['stop']   )    }
            masks.append(mask)
    
    return { 'regions' : regions, 'masks' : masks }


def write_track( data_forward, data_reverse, file, track_name ) :
    """
    Write a Gaggle Genome Browser compatible track file.
    """
    f = open( file, 'w' )
    f.write( 'sequence\tstrand\tposition\tvalue\n' )
    for n,i in enumerate( data_forward ) :
        if i != 0 :
            f.write( track_name + '\t+\t' + str(n) + '\t' + str(i) +
'\n' )
    for n,i in enumerate( data_reverse ) :
        if i != 0 :
            f.write( track_name + '\t-\t' + str(n) + '\t' + str(i) +
'\n' )
    f.close()

