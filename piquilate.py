#!/usr/bin/env python
"""
Report the amplitudes with respect to the bacground in a region list.
"""
import pique
import yaml
import sys

str_opts = [    'track_name',           \
                'annotated_bookmarks',  \
                'gene_annotations',     \
                'new_bookmarks',        ]

opt_dict = yaml.load( open( sys.argv[1] ).read() )

for opt in str_opts :
    if not opt_dict.has_key( opt ) :
        print 'config file missing option : ' + opt
        quit()
    setattr( sys.modules[__name__], opt, opt_dict[opt] )

# read bookmarks file
pique.msg( 'reading annotations...' )
peaks = pique.readbookmarks( annotated_bookmarks )

# read gene annotations
genes = {}
for line in open( gene_annotations ) :
    if line.__contains__( '\"' ) :
        continue
    if not line.split()[0].lower() == track_name :
        continue
    start,stop = map( int, line.split()[1:3] )
    strand = line.strip().split()[5]
    name = line.split()[3]
    genes[name] = { 'start':start,'stop':stop,'strand':strand }
    print name
print len(genes.keys())

# add gene annotation to peak list
pique.msg( 'finding transcript start sites in enriched regions...' )
for peak in peaks :
    gg = []
    for name,gene in genes.items() :
        if peak['start'] < gene['start'] < peak['stop'] :
            gg.append(name)
    peak['annotations']['genes'] = ','.join(gg)

# write new bookmark file
pique.msg( 'writing re-annotated bookmark file...' )
pique.writebookmarks( new_bookmarks, track_name, peaks )
