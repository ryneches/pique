#!/usr/bin/env python

import pique
import yaml
import sys

str_opts = [    'track_name',           \
                'forward_ChIP_track',   \
                'forward_bgnd_track',   \
                'reverse_ChIP_track',   \
                'reverse_bgnd_track',   \
                'slice_bookmarks',      \
                'new_track_prefix',     ]

opt_dict = yaml.load( open( sys.argv[1] ).read() )

for opt in str_opts :
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

# read bookmarks file
pique.msg( 'reading annotations...' )
slices = pique.readbookmarks( slice_bookmarks )

# write new slice tracks
for s in slices :
    sdata_ff    = data_ff[ s['start'] : s['stop'] ]
    sdata_rr    = data_rr[ s['start'] : s['stop'] ]
    sb_ff       = b_ff[    s['start'] : s['stop'] ]
    sb_rr       = b_rr[    s['start'] : s['stop'] ]
    
    # write enrichment track
    file = new_track_prefix + '_IP_' + s['annotations']['slice'] + '.track'
    pique.msg( 'writing ' + file + '...' )
    pique.write_track( sdata_ff, sdata_rr, file, track_name )

    # write background track
    file = new_track_prefix + '_BG_' + s['annotations']['slice'] + '.track'
    pique.msg( 'writing ' + file + '...' )
    pique.write_track( sb_ff, sb_rr, file, track_name )

    # zero out sliced region
    data_ff[ s['start'] : s['stop'] ] = 0
    data_rr[ s['start'] : s['stop'] ] = 0
    b_ff[    s['start'] : s['stop'] ] = 0
    b_rr[    s['start'] : s['stop'] ] = 0

# write new tracks, minus sliced regions
pique.msg( 'writing sliced tracks...' )
file = new_track_prefix + '_IP_' + track_name + '.track'
pique.write_track( data_ff, data_rr, file, track_name )
file = new_track_prefix + '_BG_' + track_name + '.track'
pique.write_track( b_ff,    b_rr,    file, track_name )

