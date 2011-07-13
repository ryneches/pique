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
        pique.msg( 'starting run for project : ' + name )
        
        alpha    = int( self.alphatext.get().strip() )
        l_thresh = int( self.lthreshtext.get().strip() )
        
        # load the data
        pique.msg( 'loading data...' )
        self.master.title( 'Pique : loading data...' )
        if not self.mapfile :
            D = pique.data.PiqueData( self.IPfile, self.BGfile )
        else :
            D = pique.data.PiqueData( self.IPfile, self.BGfile, self.mapfile )     
        
        # start analysis workbench
        pique.msg( 'creating analysis workbench...' )
        self.master.title( 'Pique : creating workbench...' )
        PA = pique.analysis.PiqueAnalysis( D )
        
        # run filters
        pique.msg( 'running filters...' )
        self.master.title( 'Pique : running filters...' )
        pique.msg( '  -> alpha    : ' + str(alpha) )        
        pique.msg( '  -> l_thresh : ' + str(l_thresh) )
        PA.filter_all( alpha, l_thresh )
        
        # find peaks
        pique.msg( 'finding peaks...' )
        self.master.title( 'Pique : finding peaks...' )
        for ar_name in PA.data.keys() :
            PA.find_peaks(ar_name)
        
        # write output files
        pique.msg( 'writing output files...' )
        self.master.title( 'Pique : writing output...' )
        pique.fileIO.writepeaksGFF( name + '.gff', PA.data )
        pique.fileIO.writebookmarks( name + '.bookmark', PA.data )
        
        # done!
        pique.msg( 'run completed.' )
        self.master.title( 'Pique : run completed.' )

app = PiqueApp( root )

root.mainloop()
