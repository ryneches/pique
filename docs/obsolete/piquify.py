#!/usr/bin/env python
import numpy
import yaml
import sys
import pique
"""
Normaize a background track relative to a ChIP track using a
bookmark file of non-peak regions.

Usage : ./piquify.py <config.yaml>
"""

# Parse config file, set runtime parameters
str_opts = [    'track_name',               \
                'forward_ChIP_track',       \
                'forward_bgnd_track',       \
                'reverse_ChIP_track',       \
                'reverse_bgnd_track',       \
                'new_forward_bgnd_track',   \
                'new_reverse_bgnd_track',   \
                'non_peak_bookmarks',       ]

opt_dict = yaml.load( open( sys.argv[1] ).read() )

for opt in str_opts :
    if not opt_dict.has_key( opt ) :
        print 'config file missing option : ' + opt
        quit()
    setattr( sys.modules[__name__], opt, opt_dict[opt] )

pique.msg( 'reading track data...' )
data_ff = pique.readtrack( forward_ChIP_track )
data_rr = pique.readtrack( reverse_ChIP_track )
b_ff    = pique.readtrack( forward_bgnd_track )
b_rr    = pique.readtrack( reverse_bgnd_track )

non_peaks = pique.readbookmarks( non_peak_bookmarks )

# calculate enrichment ratios
pique.msg( 'calculating enrichment ratios using ' + \
            str(len(non_peaks)) + 'regions...' )

d_f, d_r, b_f, b_r = [],[],[],[]
for n,region in enumerate( non_peaks ) :
    d_f.append( sum( data_ff[ region['start'] : region['stop'] ] ) )
    d_r.append( sum( data_rr[ region['start'] : region['stop'] ] ) )
    b_f.append( sum(    b_ff[ region['start'] : region['stop'] ] ) )
    b_r.append( sum(    b_rr[ region['start'] : region['stop'] ] ) )
    
f_norm = float(sum(d_f)) / sum(b_f)
r_norm = float(sum(d_r)) / sum(b_r)

pique.msg( 'forward track background scaled by ' + str(f_norm) )
pique.msg( 'reverse track background scaled by ' + str(r_norm) )

pique.msg( 'writing new background tracks...' )
f = open( new_forward_bgnd_track, 'w' )
for n,i in enumerate( b_ff * f_norm ) :
    f.write( track_name + '\t' + str(n) + '\t' + str(i) + '\n' )
f.close()

f = open( new_reverse_bgnd_track, 'w' )
for n,i in enumerate( b_rr * r_norm ) :
    f.write( track_name + '\t' + str(n) + '\t' + str(i) + '\n' )
f.close()


