#!/usr/bin/env python
import pique
import sys
import pickle
"""
Recalculates enrichment ratios from a curated bookmark file, and
prints a new bookmark file to standard output.

Usage : ./er_recalculate.py [name.pickle] [input.bookmark]
"""

PA = pickle.load( open( sys.argv[1] ) )

f = open( sys.argv[2] )
f.readline()
f.readline()

print '>name: Pique bookmarks for ' + PA.name
print 'Chromosome\tStart\tEnd\tStrand\tName\tAnnotation'

for line in f :
    contig, start, stop = line.split('\t')[:3]
    start = int(start)
    stop  = int(stop)
     
    # pick out the binding coordinate
    for annotation in line.split('\t')[5].split() :
        key,value = annotation.split(':')
        if key == 'binds_at' :
            bc = int(value)
            break
    
    # find the right analysis region (contigs may have more than one)
    test = False
    for ar in PA.data.keys() :
        if PA.data[ar]['contig'] == contig :
            if start > PA.data[ar]['region']['start'] and start < PA.data[ar]['region']['stop'] :
                test = True
                break
    if not test :
        err = 'Annotation not contained in any analysis region.\n' + line
        raise Exception( err )
    
    start = start - PA.data[ar]['region']['start']
    stop  = stop  - PA.data[ar]['region']['start']

    # calculate the new enrichment ratio
    ip_e = sum( PA.data[ar]['IP']['forward'][start:stop]    \
              + PA.data[ar]['IP']['reverse'][start:stop]    )
    bg_e = sum( PA.data[ar]['BG']['forward'][start:stop]    \
              + PA.data[ar]['BG']['reverse'][start:stop]    )
    if not float(bg_e) == 0 :
        er = float(ip_e) / float(bg_e) 
    else :
        raise Exception( 'Zero background!')
    
    annotations = []
    annotations.append( 'enrichment_ratio:' + str(er) )
    annotations.append( 'binds_at:' + str(bc) )
    for n,norm in enumerate( PA.data[ar]['norms']) :
        annotations.append( 'norm_' + str(n) + ':' + str(norm) )
    
    print '\t'.join( map( str, [ contig, start, stop, 'none', 'peak', er, ' '.join(annotations) ] ) )
    
    
