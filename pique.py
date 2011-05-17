#!/usr/bin/env python
"""
This is the main peak caller script. Once you have identified the
noise threshold in your data using piquant.py, use this script to
generate bookmark files containing putative peak annotations.

Usage : ./pique.py <config.yaml>
"""
import numpy
import sys
import yaml
import pique
import scipy.signal

# parse the config file
num_opts = [    'n_thresh', \
                'l_thresh', \
                'q_thresh', \
                'c_thresh', \
                's_thresh', \
                'too_small',\
                'too_big', ]
str_opts = [    'track_name',           \
                'forward_ChIP_track',   \
                'forward_bgnd_track',   \
                'reverse_ChIP_track',   \
                'reverse_bgnd_track',   \
                'masking_loci',         \
                'peak_bookmarks',       \
                'weed_bookmarks',       \
                'overlap_track',        \
                'binding_track' ]

opt_dict = yaml.load( open( sys.argv[1] ).read() )

for opt in num_opts + str_opts :
    if not opt_dict.has_key( opt ) :
        print 'config file missing option : ' + opt
        quit()
    setattr( sys.modules[__name__], opt, opt_dict[opt] )

# read the track data
pique.msg( 'reading track data...' )
data_ff = pique.readtrack( forward_ChIP_track )
data_rr = pique.readtrack( reverse_ChIP_track )
b_ff    = pique.readtrack( forward_bgnd_track )
b_rr    = pique.readtrack( reverse_bgnd_track )

# apply mask
pique.msg( 'applying mask...' )
is_elements = []
for line in open( masking_loci ) :
    if line.__contains__('#') :
        continue
    start, stop = map( int, line.split()[:2] )
    is_elements.append( { 'start':start, 'stop':stop } )

data_ff = pique.mask( data_ff, is_elements )
data_rr = pique.mask( data_rr, is_elements )
b_ff    = pique.mask( b_ff,    is_elements )
b_rr    = pique.mask( b_rr,    is_elements )

# run the filter
pique.msg( 'running filters...' )
data_f  = pique.filterset( data_ff, 300, l_thresh )
data_r  = pique.filterset( data_rr, 300, l_thresh )
b_f     = pique.filterset( b_ff,    300, l_thresh )
b_r     = pique.filterset( b_rr,    300, l_thresh )

# find the peaks
pique.msg( 'finding peaks...' )
forward  = pique.findregions( data_f, n_thresh )
reverse  = pique.findregions( data_r, n_thresh )
envelope = pique.overlaps( forward, reverse )

# delete peaks that aren't above the control threshold
pique.msg( 'removing peaks that are lost in the background...' )

weeds = []
for peak in envelope :
    fstart = peak['forward']['start']
    fstop  = peak['forward']['stop']
    rstart = peak['reverse']['start']
    rstop  = peak['reverse']['stop']
    a = float(sum(data_ff[fstart:fstop]) + sum(data_rr[rstart:rstop]) )
    b = float(sum(    b_f[fstart:fstop]) + sum(    b_r[rstart:rstop]) )
    if b != 0 :
        if a/b < c_thresh :
            envelope.remove(peak)
            weeds.append(peak)
            continue
    if rstop - fstart > too_big or rstop - fstart < too_small :
        envelope.remove(peak)
        weeds.append(peak)

# find putative binding coordinate in good peaks
for n,peak in enumerate( envelope ) :
    df = data_ff[ peak['forward']['start'] : peak['reverse']['stop'] ]
    dr = data_rr[ peak['forward']['start'] : peak['reverse']['stop'] ]
    df = scipy.signal.convolve( df, numpy.ones(10) / 10.0 )
    dr = scipy.signal.convolve( dr, numpy.ones(10) / 10.0 )
    c = peak['forward']['start'] + ( dr.argmax() + df.argmax() ) / 2.0
    envelope[n]['binds'] = int(c)
for n in range(len(weeds)) :
    weeds[n]['binds'] = ''

# write output files
pique.msg( 'writing output files...' )

peaks_f = numpy.zeros( len(data_f) )
peaks_r = numpy.zeros( len(data_r) )
for peak in forward :
    peaks_f[ peak['start'] : peak['stop'] ] = 1
for peak in reverse :
    peaks_r[ peak['start'] : peak['stop'] ] = 1

pique.write_strandless_track( peaks_f+peaks_r, overlap_track, track_name )

pique.writebookmarks( peak_bookmarks, track_name, envelope )
pique.writebookmarks( weed_bookmarks, track_name, weeds    )

f = open( binding_track, 'w' )
f.write( 'sequence\tstrand\tposition\tvalue\n' )
for peak in envelope :
    f.write( track_name + '\t.\t' + str(peak['binds']) + '\t1\n' )
f.close()
