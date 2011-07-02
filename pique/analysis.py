#!/usr/bin/env python
"""
Pique analysis module.
"""
import data
import processing
import numpy
import scipy


class PiqueAnalysis :
    """
    Pique analysis workbench object.
    """
    
    def __init__( self, PD ) :
        """
        Workbench initialization requires a populated PiqueData
        object. Otherwise, there's nothing to do!
        """
        self.PD = PD
        
    def noise_threshold( self, IP, BG ) :
        """
        Takes an array of IP coverage data and background coverage
        data, and searches for the noies threshold in the data.
        """
        for data in IP, BG :
            for n in arange( min(data), max(data), 
        
    def filter_data( self, contig, region, alpha, l_thresh ) :
        
        if not self.filtered.has_key( contig ) :
            self.filtered[contig] = { 'length' : self.data[contig]['length'] }
        
        for track in 'IP', 'BG' :
            
            if not self.filtered[contig].has_key( track ) :
                self.filtered[contig][track] = { 'forward' : numpy.zeros( self.data[contig]['length'] ),    \
                                                 'reverse' : numpy.zeros( self.data[contig]['length'] ) }
            
            start = region['start']
            stop  = region['stop']
            
            for strand in 'forward', 'reverse' :
                self.filtered[contig][track][strand][start:stop] =  \
                    processing.filterset( self.data[contig][track]['forward'][start:stop], alpha, l_thresh )

