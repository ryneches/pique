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
    # index the BAM file
    pysam.index( file )
    
    # open the BAM file
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
                        'stop'   : int( g['stop']  )    }
            masks.append(mask)
    
    return { 'regions' : regions, 'masks' : masks }

def writepeaksGFF( file, data ) :
    """
    Write peak list as a GFF file from a PiqueAnalysis.data object.
    """
    f = open( file, 'w' )
    for ar_name in data.keys() :
        ar      = data[ar_name]
        contig  = ar['contig']
        source  = 'Pique-1.0'
        feature = 'peak'
        strand  = '.'
        frame   = '.'
        group   = ''
        rstart  = ar['region']['start']
        for e in ar['peaks'] :
            start = str( rstart + e['start'] )
            stop  = str( rstart + e['stop']  )
            er    = str( e['annotations']['enrichment_ratio'] )
            s = '\t'.join( [contig,source,feature,start,stop,er,strand,frame,group ] )
            f.write( s + '\n' )
    f.close()
        
def writetrack( file, data, track='IP' ) :
    """
    Write a Gaggle Genome Browser compatible track file from a
    PiqueData.data object.
    """
    f = open( file, 'w' )
    f.write( 'sequence\tstrand\tposition\tvalue\n' )
    for contig in data.keys() :
        for n,i in enumerate( data[contig][track]['forward'] ) :
            if i != 0 :
                f.write( contig + '\t+\t' + str(n) + '\t' + str(i) + '\n' )
        for n,i in enumerate( data[contig][track]['reverse'] ) :
            if i != 0 :
                f.write( contig + '\t-\t' + str(n) + '\t' + str(i) + '\n' )
    f.close()

def writeQP( file, data ) :
    """
    Write GGB quantitative positional file of estimated binding
    coodinates.
    """
    f = open( file, 'w' )
    f.write( 'sequence\tstrand\tposition\tvalue\n' )
    for ar_name in data.keys() :
        ar      = data[ar_name]
        contig  = ar['contig']
        for e in ar['peaks'] :
            er    = str( e['annotations']['enrichment_ratio'] )
            bc    = str( e['annotations']['binds_at'] )
            f.write( contig + '\t'  \
                    + '.\t'         \
                    + bc + '\t'     \
                    + er + '\n'     )
    f.close()

def writebookmarks( file, data ) :
    """
    Write a GGB bookmark file from an PiqueAnalyis.data object.
    """
    f = open( file, 'w' )
    f.write( '>name: Pique bookmarks\n' )
    f.write( 'Chromosome\tStart\tEnd\tStrand\tName\tAnnotation\n' )
    for ar_name in data.keys() :
        ar      = data[ar_name]
        contig  = ar['contig']
        source  = 'Pique-1.0'
        feature = 'peak'
        strand  = '.'
        frame   = '.'
        group   = ''
        rstart  = ar['region']['start']
        for e in ar['peaks'] :
            start = str( rstart + e['start'] )
            stop  = str( rstart + e['stop']  )
            er    = str( e['annotations']['enrichment_ratio'] )
            f.write( contig + '\t'      \
                    + start + '\t'      \
                    + stop  + '\t'      \
                    + 'none\tpeak\t'    )
            if e.has_key( 'annotations' ) :
                for key,value in e['annotations'].items() :
                    f.write( str(key) + ':' + str(value) + ' ' )
            f.write( '\n' )
    f.close()


