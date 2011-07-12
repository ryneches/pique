#!/usr/bin/env python
"""
Pique.
"""
import pique
import sys

IP_file = sys.argv[1]
BG_file = sys.argv[2]

if len(sys.argv) > 3 :
    map_file = sys.argv[3]
    D = pique.data.PiqueData( IP_file, BG_file, map_file=map_file )
else :
    D = pique.data.PiqueData( IP_file, BG_file )

# FIXME : THIS IS BECAUSE OUR TEST DATA HAS ERRORS IN IT.
del(D.data['PNRC100'])

pique.msg( 'loading data...' )
PA = pique.analysis.PiqueAnalysis( D )

pique.msg( 'filtering data...' )
PA.filter_all( 300, 30 )

pique.msg( 'finding peaks...' )
for ar_name in PA.data.keys() :
    PA.find_peaks(ar_name)

pique.msg( 'writing output files...' )
pique.fileIO.writepeaksGFF( 'test_out.gff', PA.data )
pique.fileIO.writebookmarks( 'test.bookmark', PA.data )

