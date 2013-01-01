#!/usr/bin/env python

# Pique version
version         = '0.1.5'

# cache directory for storing temporary and downloaded files
cache_dir       = '/tmp/pique'

# Test data, IP BAM file
test_ip_file    = 'http://files.figshare.com/224619/IP.bam'

# Test data, background BAM file
test_bg_file    = 'http://files.figshare.com/224618/WCE.bam'

# Test data, mapping file
test_map_file   = 'http://files.figshare.com/224617/map.gff'

# Permitted separation between peaks found on the forward and reverse
# strand. This should correspond with the difference between the
# fragment length and the read length.
peak_separation = 120

# The fall-off from maximum needed to indentify a local maximum
top_delta       = 0.01

# The fall-off from an identified local maximum used to identify the
# peak region.
reg_delta       = 0.0625

# The maximum radius around an identified local maximum that can be
# included in a peak region.
radius          = 300
