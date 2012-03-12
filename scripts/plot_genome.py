#!/usr/bin/env python
import pique
import pylab

print( 'loading data...' )
D = pique.data.PiqueData(                                                           \
    '/home/russell/Dropbox/Peaks/tfbD1_data/final.I.3171780.sam.sorted.bam',        \
    '/home/russell/Dropbox/Peaks/tfbD1_data/ctrl.final.I.3171780.sam.sorted.bam',   \
    '/home/russell/Dropbox/Peaks/tfbD1_data/new2.Hb_NRC1_mapping.gff',              \
    name='tfbD1' )

print( '  found contigs :' )
for contig in D.data.keys() :
    print( '    ' + contig )
    print( '      length : ' +
    str(D.data[contig]['length']) )
    for r in D.data[contig]['regions'] :
        start = str( r['start'] )
        stop  = str( r['stop']  )
        print( '      analysis region : ' + start + ':' + stop )
    for m in D.data[contig]['masks'] :
        start = str( m['start'] )
        stop  = str( m['stop']  )
        print( '      masking region  : ' + start + ':' + stop )
    
# start analysis workbench
PA = pique.analysis.PiqueAnalysis( D )

print( '   applying filters...' )
for ar_name in PA.data.keys() :
    PA.apply_filter( ar_name, 30, 300 )

# find peaks
print( '   finding peaks...' )
for ar_name in PA.data.keys() :
    PA.find_peaks(ar_name)

# load CSDeconv regions
f = open( '/home/russell/Projects/pique/docs/Csdeconv-TfbD-ER.txt' )
f.readline()
csd_regions = []
for line in f :
    contig, start, stop, a_1, a_2 = line.split()
    if contig == 'Chromosome' :
        start, stop = int(start), int(stop)
        csd_regions.append( { 'start' : start, 'stop' : stop } )
f.close()

# load transcript start sites and predicted binding motifs
tsses = []
motifs = []
for line in open( '/home/russell/Projects/pique/docs/ExpTSSandMotifs.txt' ) :
    contig,tss,mstart,mstop = line.split()
    if contig.split(':')[0] == 'Chromosome' :
        if not tss == '.' :
            tsses.append( int(tss) )
        if not mstart == '.' :
            motifs.append( { 'start' : int(mstart), 'stop' : int(mstop) } )
            mstart,mstop = int(mstart), int(mstop)

# mark peak regions
for peak in PA.data['Chromosome_0:2014239']['peaks'] :
    axvspan( peak['start'], peak['stop'], color='orange', alpha=0.3 )

for peak in csd_regions :
    axvspan( peak['start'], peak['stop'], color='blue', alpha=0.3 )

# mark transcript start sites
for tss in tsses :
    pylab.arrow( tss, -10, 0, 0, width=5, color='black' )

# mark predicted binding motifs
for motif in motifs :
    m = motif['start'] + ( motif['stop'] - motif['start'] ) / 2.0
    pylab.arrow( m, -10, 0, 0, width=5, color='cyan' )

# place gene annotations
genes = {}
for line in open( '/home/russell/Dropbox/Peaks/genes/annotations.bed' ) :
    if line.__contains__('Chromosome') :
        start, stop = map( int, line.split()[1:3] )
        name = line.split()[3]
        strand = line.split()[5]
        genes[name] = {'start':start,'stop':stop,'strand':strand}

for name,region in genes.items() :
    start = region['start']
    stop  = region['stop']
    l     = stop - start
    if region['strand'] == '-' :
        #pylab.arrow( start, -10, l, 0, width=5, color='green' )
        bbox_props = dict(boxstyle="rarrow,pad=0.3", fc="green", ec="b", lw=0, alpha=0.6)
        text( start + (stop-start)/2.0, -10, name, ha="center", va="center", size=15, bbox=bbox_props)
    else :
        #pylab.arrow( start, -20, l, 0, width=5, color='purple' )
        bbox_props = dict(boxstyle="larrow,pad=0.3", fc="purple", ec="b", lw=0, alpha=0.6)
        text( start + (stop-start)/2.0, -20, name, ha="center", va="center", size=15, bbox=bbox_props)

# plot data
PA.data['Chromosome_0:2014239']['IP']['forward'][0] = 0.0 # start fill off at zero
PA.data['Chromosome_0:2014239']['IP']['reverse'][0] = 0.0 # start fill off at zero
fill( PA.data['Chromosome_0:2014239']['IP']['forward'], 'blue', alpha=0.6, lw=0 )
fill( PA.data['Chromosome_0:2014239']['IP']['reverse'], 'red', alpha=0.6, lw=0 )

plot( PA.data['Chromosome_0:2014239']['ip']['forward'], 'blue', lw=2 )
plot( PA.data['Chromosome_0:2014239']['ip']['reverse'], 'red', lw=2 )

# plot threshold
axhline( PA.data['Chromosome_0:2014239']['n_thresh'], color='black', ls='--', lw=2 )
