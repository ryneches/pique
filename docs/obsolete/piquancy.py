#!/usr/bin/env python
"""
Report the amplitudes with respect to the bacground in a region list.
"""
import pique
import yaml
import sys

str_opts = [    'track_name',           \
                'peak_bookmarks',       \
                'forward_ChIP_track',   \
                'forward_bgnd_track',   \
                'reverse_ChIP_track',   \
                'reverse_bgnd_track',   \
                'masking_loci',         \
                'annotated_bookmarks'   ]

opt_dict = yaml.load( open( sys.argv[1] ).read() )

for opt in str_opts :
    if not opt_dict.has_key( opt ) :
        print 'config file missing option : ' + opt
        quit()
    setattr( sys.modules[__name__], opt, opt_dict[opt] )

# read track data
pique.msg( 'reading track data...' )
data_ff = pique.readtrack( forward_ChIP_track )
data_rr = pique.readtrack( reverse_ChIP_track )
b_ff    = pique.readtrack( forward_bgnd_track )
b_rr    = pique.readtrack( reverse_bgnd_track )

# read bookmarks file
peaks = pique.readbookmarks( peak_bookmarks )

# calculate enrichment ratios
for n,peak in enumerate(peaks) :
    a =     sum( data_ff[ peak['start'] : peak['stop'] ] )
    a = a + sum( data_rr[ peak['start'] : peak['stop'] ] )
    b =     sum(    b_ff[ peak['start'] : peak['stop'] ] )
    b = b + sum(    b_rr[ peak['start'] : peak['stop'] ] )
    peaks[n]['annotations']['enrichment_ratio'] = float(a) / float(b)

# write new bookmark file
pique.msg( 'writing re-annotated bookmark file...' )
pique.writebookmarks( annotated_bookmarks, track_name, peaks )
