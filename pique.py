#!/usr/bin/env python
"""
Pique.
"""
import pique
import sys
from Tkinter import *
from tkFileDialog import *

# create the main window
root = Tk()

class PiqueApp :
    def __init__( self, master ) :
        self.IPfile  = None
        self.BGfile  = None
        self.mapfile = None
        self.alpha   = 300
        self.l_thresh= 30
        self.name    = 'MyProject'
        self.master = master

        master.title( 'Pique 1.0' )
        
        # project name
        frame0 = Frame( master )
        self.nametext     = Entry(  frame0                  )
        self.nametext.pack( side=LEFT )
        self.nametext.insert( 0, self.name )
        
        self.namelabel = Label( frame0, text='Project name' )
        self.namelabel.pack( side=RIGHT, fill=X )
        frame0.pack( fill=X )

        # IP file 
        frame1 = Frame( master )
        self.IPfilebutton = Button( frame1,                 \
                                    text='IP Track',        \
                                    command=self.setIPfile  )
        self.IPfilebutton.pack( side=RIGHT, fill=X )
        
        self.IPfiletext  = Entry(   frame1                  )
        self.IPfiletext.pack( side=LEFT )
        frame1.pack( fill=X )
         
        # BG file 
        frame2 = Frame( master )
        self.BGfilebutton = Button( frame2,                 \
                                    text='BG Track',        \
                                    command=self.setBGfile  )
        self.BGfilebutton.pack( side=RIGHT, fill=X )
        
        self.BGfiletext  = Entry(   frame2                  )
        self.BGfiletext.pack( side=LEFT )
        frame2.pack( fill=X )
        
        # map file
        frame3 = Frame( master )
        self.mapfilebutton = Button(frame3,                 \
                                    text='Map (optional)',  \
                                    command=self.setmapfile )
        self.mapfilebutton.pack( side=RIGHT, fill=X )
        
        self.mapfiletext  = Entry(  frame3                  )
        self.mapfiletext.pack( side=LEFT )
        frame3.pack( fill=X )
        
        # alpha
        frame4 = Frame( master )
        self.alphatext    = Entry(  frame4                  )
        self.alphatext.pack( side=LEFT )
        self.alphatext.insert( 0, self.alpha )        
        
        self.alphalabel = Label( frame4, text='Fragment length' )
        self.alphalabel.pack( side=RIGHT )
        frame4.pack( fill=X )       
        
        # l_thresh
        frame5 = Frame( master )
        self.lthreshtext  = Entry(  frame5                  )
        self.lthreshtext.pack( side=LEFT )
        self.lthreshtext.insert( 0, self.l_thresh )        
        
        self.lthreshlabel = Label( frame5, text='Read length' )
        self.lthreshlabel.pack( side=RIGHT )
        frame5.pack( fill=X )
        
        # control buttons
        frame6 = Frame( master )
        self.runbutton    = Button( frame6,                 \
                                    text='Run',             \
                                    command=self.run        )
        self.runbutton.pack( side=RIGHT )
        
        self.quitbutton   = Button( frame6,                 \
                                    text='Quit',            \
                                    command=master.quit     )
        self.quitbutton.pack( side=LEFT )
        frame6.pack( fill=X )
        
    def setIPfile( self ) :
        self.IPfile = askopenfilename()
        l = len(self.IPfiletext.get())
        self.IPfiletext.delete( 0, l )
        self.IPfiletext.insert( 0, self.IPfile )
        print 'IP file : ' + self.IPfile
        
    def setBGfile( self ) :
        self.BGfile = askopenfilename()
        l = len(self.BGfiletext.get())
        self.BGfiletext.delete( 0, l )
        self.BGfiletext.insert( 0, self.BGfile )
        print 'BG file : ' + self.BGfile
        
    def setmapfile( self ) :
        self.mapfile = askopenfilename()
        l = len(self.mapfiletext.get())
        self.mapfiletext.delete( 0, l )
        self.mapfiletext.insert( 0, self.mapfile )
        print 'map file : ' + self.mapfile
        
    def run( self ) :
        
        # check inputs...
        name     = self.nametext.get().strip()
        
        # set logfile
        logfile = name + '.log'
        
        pique.msg( logfile, 'starting run for project : ' + name )
         
        alpha    = int( self.alphatext.get().strip() )
        l_thresh = int( self.lthreshtext.get().strip() )
               
        # log inputs
        pique.msg( logfile, '  -> IP file  : ' + self.IPfile   )
        pique.msg( logfile, '  -> BG file  : ' + self.BGfile   )
        pique.msg( logfile, '  -> map file : ' + self.mapfile  )
        pique.msg( logfile, '  -> alpha    : ' + str(alpha)    )
        pique.msg( logfile, '  -> l_thresh : ' + str(l_thresh) )
        
        # load the data
        pique.msg( logfile, 'loading data...' )
        self.master.title( 'Pique : loading data...' )
        if not self.mapfile :
            D = pique.data.PiqueData( self.IPfile, self.BGfile, name=name )
        else :
            D = pique.data.PiqueData( self.IPfile, self.BGfile, self.mapfile, name=name )
        
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
        self.master.title( 'Pique : creating workbench...' )
        PA = pique.analysis.PiqueAnalysis( D )
        
        # run filters
        pique.msg( logfile, 'running filters...' )
        self.master.title( 'Pique : running filters...' )
        
        for ar_name in PA.data.keys() :
            pique.msg( logfile, '  :: applying filters to analysis region ' + ar_name )
            PA.apply_filter( ar_name, alpha, l_thresh )
            
        # find peaks
        pique.msg( logfile, 'finding peaks...' )
        self.master.title( 'Pique : finding peaks...' )
        for ar_name in PA.data.keys() :
            PA.find_peaks(ar_name)
            pique.msg( logfile, '  peaks ' + ar_name + ' : ' + str(len(PA.data[ar_name]['peaks'])) )
            pique.msg( logfile, '     noise threshold : ' + str(PA.data[ar_name]['N_thresh']) )
            pique.msg( logfile, '     normalizations  : ' + ', '.join( map(str, PA.data[ar_name]['norms']) ) )

        # write output files
        pique.msg( logfile, 'writing output files...' )
        self.master.title( 'Pique : writing output...' )
        pique.fileIO.writepeaksGFF(  name + '.gff',      PA.data )
        pique.fileIO.writebookmarks( name + '.bookmark', PA.data, name=name )
        pique.fileIO.writeQP(        name + '.qp',       PA.data )
        pique.fileIO.writepeakTSV(   name + '.peak.tsv', PA.data )
        pique.fileIO.writetrack(     name + '.IP.track', D.data  )
        pique.fileIO.writetrack(     name + '.BG.track', D.data, track='BG' )
        
        # done!
        pique.msg( logfile, 'run completed.' )
        self.master.title( 'Pique : run completed.' )

app = PiqueApp( root )

root.mainloop()
