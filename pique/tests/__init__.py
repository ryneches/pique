"""
nose tests for pique.
"""
import pique.data
import pique.constants
import numpy

def test_init_PiqueData() :
    
    # does the data load?
    D = pique.data.PiqueData( pique.constants.test_ip_file,     \
                              pique.constants.test_bg_file,     \
                              pique.constants.test_map_file,    \
                              name='test' )
    
    # are they the ones we were looking for?
    assert D.data.keys().__contains__( 'Chromosome' )
    assert D.data.keys().__contains__( 'PNRC200' )
    
    # is the track data present?
    for contig in D.data.keys() :
        assert type(D.data[contig]['IP']['forward']) is numpy.ndarray
        assert type(D.data[contig]['IP']['reverse']) is numpy.ndarray
        assert type(D.data[contig]['BG']['forward']) is numpy.ndarray
        assert type(D.data[contig]['BG']['reverse']) is numpy.ndarray

def test_import_cython_stuff() :
    import pique.peak
