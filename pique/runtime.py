#!/usr/bin/env python
"""
Pique runtime.
"""
import pique
import sys
import cPickle
import numpy

def makemap( name, 
             bamfile,
             window,
             stride, 
             highest,
             lowest,
             bins ) :
    """
    This function drives the genome map making workflow.
    """
    import pylab 

    logfile = name + '.mapmaker.log'
    mapfile = name + '.map.gff'
    
    pique.msg( logfile, 'starting mapmaker for project : ' + name )
    
    pique.msg( logfile, '  -> BAM file    : ' + bamfile      )
    pique.msg( logfile, '  -> map file    : ' + mapfile      )
    pique.msg( logfile, '  -> window      : ' + str(window)  )
    pique.msg( logfile, '  -> stride      : ' + str(stride)  )
    pique.msg( logfile, '  -> bins        : ' + str(bins)    )
    pique.msg( logfile, '  -> highest bin : ' + str(highest) )
    pique.msg( logfile, '  -> lowest bin  : ' + str(lowest)  )
    
    pique.msg( logfile, 'loading data...' )
    
    data = pique.fileIO.loadBAM( bamfile )
    
    pique.msg( logfile, '  found contigs :' )
    
    for contig in data.keys() :
        pique.msg( logfile, '    ' + contig )
        pique.msg( logfile, '      ' + str(len(data[contig]['forward'])) )
    
    pique.msg( logfile, 'making spectral histograms...' )
    
    sh = {}
    for contig in data.keys() :
        pique.msg( logfile, '  :: making sectral histogram for contig ' + contig )
        d = numpy.array( data[contig]['forward'] + data[contig]['reverse'], dtype = int )
        sh[contig] = pique.mapmaker.hist(   d,
                                            lowest,
                                            highest,
                                            bins,
                                            window,
                                            stride     )
    
    # save images of spectral histograms
    pique.msg( logfile, 'saving images of spectral histograms...' )
    
    for contig in sh.keys() :
        pylab.cla() # clean up crumbs from last plot
        pylab.clf() # clean up crumbs from last plot
        pique.msg( logfile, '  :: saving image for contig ' + contig )
        pylab.contourf( sh[contig], bins )
        pylab.title( name + ' : ' + contig )
        imgname = name + '_' + contig + '.png'
        pylab.savefig( imgname, format='png' )
    

def detect( name,
            ipfile, 
            bgfile,
            mapfile,
            alpha,
            l_thresh,
            pickle_file,
            wav_file ) :
    """
    This function drives the peak detection workflow.
    """
    # set logfile
    logfile = name + '.log'
    
    pique.msg( logfile, 'starting run for project : ' + name )
    
    # log inputs
    pique.msg( logfile, '  -> IP file  : ' + ipfile   )
    pique.msg( logfile, '  -> BG file  : ' + bgfile   )
    pique.msg( logfile, '  -> map file : ' + mapfile  )
    pique.msg( logfile, '  -> alpha    : ' + str(alpha)    )
    pique.msg( logfile, '  -> l_thresh : ' + str(l_thresh) )
    
    # load the data
    pique.msg( logfile, 'loading data...' )
    D = pique.data.PiqueData( ipfile, bgfile, mapfile, name=name )
    
    pique.msg( logfile, '  found contigs :' )
    for contig in D.data.keys() :
        pique.msg( logfile, '    ' + contig )
        pique.msg( logfile, '      length : ' + str(D.data[contig]['length']) )
        for r in D.data[contig]['regions'] :
            start = str( r['start'] )
            stop  = str( r['stop']  )
            pique.msg( logfile, '      analysis region : ' + start + ':' + stop )
        for m in D.data[contig]['masks'] :
            start = str( m['start'] )
            stop  = str( m['stop']  )
            pique.msg( logfile, '      masking region  : ' + start + ':' + stop )
    
    # start analysis workbench
    pique.msg( logfile, 'creating analysis workbench...' )
    PA = pique.analysis.PiqueAnalysis( D )
    
    # run filters
    pique.msg( logfile, 'running filters...' )
    
    for ar_name in PA.data.keys() :
        pique.msg( logfile, '  :: applying filters to analysis region ' + ar_name )
        PA.apply_filter( ar_name, alpha, l_thresh )
    
    # find peaks
    pique.msg( logfile, 'finding peaks...' )
    for ar_name in PA.data.keys() :
        PA.find_peaks(ar_name)
        pique.msg( logfile, '  peaks ' + ar_name + ' : ' + str(len(PA.data[ar_name]['peaks'])) )
        pique.msg( logfile, '     noise threshold  : ' + str(PA.data[ar_name]['N_thresh']) )
        pique.msg( logfile, '     filter threshold : ' + str(PA.data[ar_name]['n_thresh']) )
        pique.msg( logfile, '     normalizations   : ' + ', '.join( map(str, PA.data[ar_name]['norms']) ) )
    
    # if a pickle file was requested, write it
    if pickle_file :
        pique.msg( logfile, 'pickling analysis workbench...' )
        cPickle.dump( PA, open( name + '.pickle', 'w' ) )
    
    # if a WAV file was requested, write it
    if wav_file :
        for contig in D.data.keys() :
            file = name + '_' + contig + '.wav'
            pique.msg( logfile, 'writing WAV output : ' + file )
            pique.fileIO.writeWAV(  file,
                                    D.data,
                                    contig,
                                    track='IP',
                                    minusBG=True,
                                    amplify=True )
    
    
    # write output files
    pique.msg( logfile, 'writing output files...' )
    pique.fileIO.writepeaksGFF(  name + '.gff',      PA.data )
    pique.fileIO.writebookmarks( name + '.bookmark', PA.data, name=name )
    pique.fileIO.writeQP(        name + '.qp',       PA.data )
    pique.fileIO.writepeakTSV(   name + '.peak.tsv', PA.data )
    pique.fileIO.writetrack(     name + '.IP.track', D.data  )
    pique.fileIO.writetrack(     name + '.BG.track', D.data, track='BG' )

    # done!
    pique.msg( logfile, 'run completed.' )

def bam2wav( name,
             ipfile, 
             bgfile ) :
    """
    This function drives the creation of a WAV file from a BAM file.
    """
    
    # set logfile
    logfile = name + '.log'
    
    pique.msg( logfile, 'converting BAM files to WAV files for project : ' + name )
    
    # log inputs
    pique.msg( logfile, '  -> IP file  : ' + ipfile   )
    pique.msg( logfile, '  -> BG file  : ' + bgfile   )
    
    # load the data
    pique.msg( logfile, 'loading data...' )
    D = pique.data.PiqueData( ipfile, bgfile, '', name=name )
    
    pique.msg( logfile, '  found contigs :' )
    for contig in D.data.keys() :
        pique.msg( logfile, '    ' + contig )
        pique.msg( logfile, '      length : ' + str(D.data[contig]['length']) )
        for r in D.data[contig]['regions'] :
            start = str( r['start'] )
            stop  = str( r['stop']  )
            pique.msg( logfile, '      analysis region : ' + start + ':' + stop )
        for m in D.data[contig]['masks'] :
            start = str( m['start'] )
            stop  = str( m['stop']  )
            pique.msg( logfile, '      masking region  : ' + start + ':' + stop )
    # write the WAV files
    for contig in D.data.keys() :
        file = name + '_' + contig + '.wav'
        pique.msg( logfile, 'writing WAV output : ' + file )
        pique.fileIO.writeWAV( file,
                               D.data,
                               contig,
                               track='IP',
                               minusBG=True,
                               amplify=True )

    # done!
    pique.msg( logfile, 'conversion completed.' )


