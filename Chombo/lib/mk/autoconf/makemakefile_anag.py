#
# Run by bootstrap script.  PWD is always the top src dir (i.e. where bootstrap is).
#

import sys
import os

subdir=sys.argv[1]
subdirdepth=subdir.count('/')

path_to_anag_style_build = subdirdepth * '../' + 'lib/mk/autoconf/anag_style_build.py'
python_cmd = '\t' + 'python ' + path_to_anag_style_build + ' ' \
                + '$(MAKE) "|||" $(MAKEFLAGS)'

sys.stdout.write( '.PHONY: install all clean uninstall realclean\n' )

for target in ( 'install', 'all', 'clean', 'uninstall', 'realclean') :
    sys.stdout.write( target + ':\n' )
    sys.stdout.write( python_cmd + ' ' + target + '\n')
