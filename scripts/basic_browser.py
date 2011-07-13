#!/usr/bin/env python
"""
A very, very basic genome browser.

usage : ./basic_browser.py forward.track    \
                           reverse.track    \
                           peak.bookmark    \
                           genes.bed
"""
import numpy
import pique
import pylab
import sys

# read track data
pique.msg( 'reading track data...' )
data_ff = pique.readtrack( sys.argv[1] )
data_rr = pique.readtrack( sys.argv[2] )

# read peak data
pique.msg( 'reading peak bookmarks...' )
peaks = pique.readbookmarks( sys.argv[3] )

# draw peaks 
pique.msg( 'drawing peak bookmarks...' )
for peak in peaks :
    axvspan( peak['start'], peak['stop'], color='green', alpha=0.3 )

# read BED formatted gene annotations
pique.msg( 'reading gene annotations...' )
genes = {}
for line in open( sys.argv[4] ) :
    if line.__contains__('Chromosome') :
        start, stop = map( int, line.split()[1:3] )
        name = line.split()[3]
        strand = line.split()[5]
        genes[name] = {'start':start,'stop':stop,'strand':strand}

# draw gene annotations
pique.msg( 'drawing gene annotations...' )
for name,region in genes.items() :
    start = region['start']
    stop  = region['stop']
    xbox  = numpy.array( [ start,   start,  stop,   stop    ] )
    ybox  = numpy.array( [ 0,       50,     50,     0       ] )
    if region['strand'] == '-' :
        ybox = -ybox
        fill( xbox, ybox, 'r', alpha=0.5 )
        text( start, -60, name )
    else :
        fill( xbox, ybox, 'b', alpha=0.5 )
        text( start, 60, name)

show()
