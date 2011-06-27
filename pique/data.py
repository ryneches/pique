#!/usr/bin/env python
"""
Data representations
"""
import numpy
import fileIO

class PiqueDataException( Exception ) :
    pass

class PiqueData :
    """
    Container class for managing Pique data.
    """
    def __init__( self ) :
        self.data = {}

    def add_contig( self,   contig_name,    \
                            IP_forward,     \
                            IP_reverse,     \
                            BG_forward,     \
                            BG_reverse      ) :
        
        l = map(len, [IP_forward, IP_reverse, BG_forward, BG_reverse ] )
        if not all( x == l[0] for x in l ) :
            raise PiqueDataException( 'Track have different lengths.' )

        IP = { 'forward' : IP_forward, 'reverse' : IP_reverse }
        BG = { 'forward' : BG_forward, 'reverse' : BG_reverse }

        self.data[ contig_name ] = { 'IP' : IP, 'BG' : BG }
        
        self.data[ contig_name ][ 'length' ] = len( IP_forward )
    
    def load_data( self, IP_file, BG_file ) :
        """
        Digest and load data from BAM files containing alignments of
        IP and background reads.

        BAM files must be prepared in such a way that the contig names
        are identical, and each IP track must be identical in length
        to its corresponding background contig.
        """
        
        IP_tracks = fileIO.loadBAM( IP_file )
        BG_tracks = fileIO.loadBAM( BG_file )
        
        IP_contigs = IP_tracks.keys().sort()
        BG_contigs = BG_tracks.keys().sort()
        
        if not len(IP_contigs) == len(BG_contigs) :
            raise PiqueDataException( 'BG and IP have different number of contigs.' )
        
        if not all( map( lambda x : x[0] == x[1], zip(IP_contigs,BG_contigs) ) ) :
            raise PiqueDataException( 'BG and IP contig names do not match.',   \
                                    { 'IP' : IP_contigs, 'BG' : BG_contigs } )
        
        for contig in IP_contigs :
            IP_forward = IP_trakcs[contig]['forward']
            IP_reverse = IP_tracks[contig]['reverse']
            BG_forward = BG_tracks[contig]['forward']
            BG_reverse = BG_tracks[contig]['reverse']
            
            self.add_contig( contig,    IP_forward, \
                                        IP_reverse, \
                                        BG_forward, \
                                        BG_reverse )


