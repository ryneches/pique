#!/usr/bin/env python
"""
Pique analysis module.
"""
import data
import processing
import numpy
import scipy
import pique
import stats


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
        
        self.data = {}
        
        # populate local data container
        
        for contig in self.PD.data.keys() :
            
            ip_f = self.PD.data[contig]['IP']['forward']
            ip_r = self.PD.data[contig]['IP']['reverse']
            bg_f = self.PD.data[contig]['BG']['forward']
            bg_r = self.PD.data[contig]['BG']['reverse']
                    
            for mask in self.PD.data[contig]['masks'] :
                l = mask['stop'] - mask['start']
                
                ip_f[ mask['start'] : mask['stop'] ] = numpy.zeros(l)
                ip_r[ mask['start'] : mask['stop'] ] = numpy.zeros(l)
                bg_f[ mask['start'] : mask['stop'] ] = numpy.zeros(l)
                bg_r[ mask['start'] : mask['stop'] ] = numpy.zeros(l)
            
            for r in self.PD.data[contig]['regions'] :
                
                IP = { 'forward' : ip_f[ r['start'] : r['stop'] ],  \
                       'reverse' : ip_r[ r['start'] : r['stop'] ]   } 
                BG = { 'forward' : bg_f[ r['start'] : r['stop'] ],  \
                       'reverse' : bg_r[ r['start'] : r['stop'] ]   } 
                
                ar = { 'contig' : contig, 'IP' : IP, 'BG' : BG, 'region' : r }
                
                ar['N_thresh'] = self.noise_threshold( ar['BG']['forward'] + ar['BG']['reverse'] )

                name = contig + '_' + str( r['start'] ) + ':' + str( r['stop'] )
                
                self.data[name] = ar
                
    def noise_threshold( self, data ) :
        """
        Computes the noise threshold in an analysis region. For now,
        this is the 90th quantile of the data.
        """
        return stats.scoreatpercentile( data.tolist(), 90 )
    
    def apply_filter( self, ar_name, alpha, l_thresh, ) :
        """
        Apply the filter set to an analysis region.
        """
        ar = self.data[ar_name]

        fip_f = processing.filterset( ar['IP']['forward'], alpha, l_thresh )
        fip_r = processing.filterset( ar['IP']['reverse'], alpha, l_thresh )
        fbg_f = processing.filterset( ar['BG']['forward'], alpha, l_thresh )
        fbg_r = processing.filterset( ar['BG']['reverse'], alpha, l_thresh )
        
        self.data[ar_name]['ip'] = { 'forward' : fip_f, 'reverse' : fip_r }
        self.data[ar_name]['bg'] = { 'forward' : fbg_f, 'reverse' : fbg_r }
        self.data[ar_name]['n_thresh'] = self.noise_threshold( fbg_f, fbg_r )
    
    def filter_all( self, alpha, l_thresh ) :
        for ar_name in self.data.keys() :
            pique.msg( ':: applying filters to analysis region ' + ar_name )
            self.apply_filter( ar_name, alpha, l_thresh )

