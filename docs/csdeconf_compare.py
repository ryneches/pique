#!/usr/bin/env python
import pylab

f = open( 'tfbd1.peak.tsv' )
f.readline()

contig,start,stop,binds_at,er = [],[],[],[],[]
for line in f :
    a,b,c,d,e = line.split()[:5]
    contig.append(      a )
    start.append(       int(b) )
    stop.append(        int(c) )
    binds_at.append(    float(d) )
    er.append(          float(e) )

f.close()
