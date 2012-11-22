#!/usr/bin/env python

# Permitted separation between peaks found on the forward and reverse
# strand. This should correspond with the difference between the
# fragment length and the read length.
peak_separation = 50

# The fall-off from maximum needed to indentify a local maximum
top_delta       = 0.05

# The fall-off from an identified local maximum used to identify the
# peak region.
reg_delta       = 0.125

# The maximum radius around an identified local maximum that can be
# included in a peak region.
radius          = 1000
