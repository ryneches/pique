#!/usr/bin/env python2.6
"""
"""
import scipy.signal
import numpy
import itertools
import operator
import sys
import warnings

# suppress annoying scypy warnings
warnings.simplefilter("ignore", DeprecationWarning)

class PiqueException( Exception ) :
    pass

def readtrack( filename ) :
    """
    Read a track file, return a numpy array.
    """
    data = []
    for line in open( filename ) :
        if line.__contains__( '\"' ) :
            continue
        data.append( float( line.split()[2] ) )
    if len( data ) == 0 :
        raise PiqueException( 'Empty input file or parsing error : ' + filename )
    return numpy.array( data )

def write_strandless_track( data, file, track_name ) :
    f = open( file, 'w' )
    f.write( 'sequence\tstrand\tposition\tvalue\n' )
    for n,i in enumerate( data ) :
        if i != 0 :
            f.write( track_name + '\t.\t' + str(n) + '\t' + str(i) + '\n' )
    f.close()

def write_track( data_forward, data_reverse, file, track_name ) :
    f = open( file, 'w' )
    f.write( 'sequence\tstrand\tposition\tvalue\n' )
    for n,i in enumerate( data_forward ) :
        if i != 0 :
            f.write( track_name + '\t+\t' + str(n) + '\t' + str(i) + '\n' )
    for n,i in enumerate( data_reverse ) :
        if i != 0 :
            f.write( track_name + '\t-\t' + str(n) + '\t' + str(i) + '\n' )
    f.close()

def read_track( file, track_name ) :
    data_ff = []
    data_rr = []
    for line in open( filename ) :
        if line.__contains__( '\"' ) :
            continue
        if not line.__contains__( track_name ) :
            continue
        strand = line.split()[1]
        if strand == '+' :
            data_ff.append( float( line.split()[3] ) )
        if strand == '-' :
            data_rr.append( float( line.split()[3] ) )
    return {'forward':numpy.array(data_ff), 'reverse':numpy.array(data_rr) }

def readbookmarks( filename ) :
    """
    Read a GGB bookmark file.
    """
    regions = []
    for line in open( filename ) :
        if line.__contains__('>') :
            continue
        if line.__contains__('\tName\tAnnotation') :
            continue
        chromosome,start,stop,strand,name,annotations = line.strip().split('\t')
        start, stop = map( int, (start,stop) )
        annot = {}
        for annotation in annotations.split() :
            if annotation.__contains__(':') :
                key,value = annotation.split(':')
                annot[key] = value
        regions.append( {   'start':start,          \
                            'stop':stop,            \
                            'strand':strand,        \
                            'chromosome':chromosome,\
                            'annotations':annot     } )
    return regions

def writebookmarks( filename, track_name, regions ) :
    """
    Write a GGB bookmark file, where any keys besides 'forward' and
    'reverse' are added to the annotation column.
    """
    f = open( filename, 'w' )
    f.write( '>name: ' + track_name + ' peak bookmarks\n' )
    f.write( 'Chromosome\tStart\tEnd\tStrand\tName\tAnnotation\n' )
    for r in regions :
        f.write(    track_name + '\t'        \
                    + str(r['start']) + '\t' \
                    + str(r['stop'] ) + '\t' \
                    + 'none\tpeak\t'         )
        if r.has_key( 'annotations' ) :
            for key,value in r['annotations'].items() :
                f.write( str(key) + ':' + str(value) + ' ' )
        f.write( '\n' )
    f.close()

def mask( data, mask_regions ) :
    """
    Mask regions in an array.
    """
    for region in mask_regions :
        data[ region['start'] : region['stop'] ] = 0
    return data

def sizeselect( regions, too_big, too_small ) :
    for region in regions :
        size = region['stop'] - region['start']
        if size > too_big or size < too_small :
            regions.remove( region )
    return regions

def msg( message ) :
    sys.stderr.write( message + '\n' )

