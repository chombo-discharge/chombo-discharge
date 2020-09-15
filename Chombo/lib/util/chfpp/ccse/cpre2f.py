# Converts ccse-style cpre files to .f files.  In the ccse trunk, this
# is done with the strip72 Perl script.  The difference here is that we
# want to skip compilation of any *_ND.F files where N != CH_SPACEDIM.
# So in those cases, we generate an empty .f file.

import os
import sys
import string

cpre     = sys.argv[1]
outfile  = sys.argv[2]
spacedim = sys.argv[3]

stump = cpre.split('.')[0]
if( (stump[-1] == 'D')
and (stump[-2] in string.digits)
and (stump[-2] != spacedim) ):
    os.system( 'echo "" > ' + outfile )
else:
    os.system( 'cat ' + cpre + ' | perl '
              + os.path.dirname(sys.argv[0]) + '/strip72 -c > ' + outfile )
