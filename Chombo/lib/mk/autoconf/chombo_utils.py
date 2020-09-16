import sys
import re

def echo_n( str ):
    """ Like 'echo -n' """
    sys.stdout.write( str )


# Arg src_test should be either 'src' or 'test', to indicate whether we are
# generating Makefile.am's under lib/src or under lib/test; the requirements
# differ a little.
def getAM_LDFLAGS( curdir, src_test ):
    """
    Print out the appropriate AM_LDFLAGS for this directory, determining it
    from the information in amldflags.dat.
    Usage: getAM_LDFLAGS('EBTools')
    """
    table = {}
    throwaway_line = re.compile('^ *\t*#|^ *\t* $|^$')

    if   src_test == 'src':
        lines = open( '../amldflags.dat' ).readlines()
    elif src_test == 'test':
        lines = open( '../amldflags.dat' ).readlines()
    else:
        assert( 0 )
    for line in lines:
        if not throwaway_line.search( line ):
            key = line.split()[0]
            vals = line.split()[1:]
            literal_vals = []
            for val in vals:
                if val in table.keys():
                    literal_vals += table[val]
                else:
                    literal_vals.append( val )
            table[key] = literal_vals
            #print "literal_vals=", literal_vals

    assert( curdir in table.keys() )

    # Now remove duplicates, add Makefile.am-ish suffixes and remove
    # Python punctuation.
    outstr = ''
    uniques = {}
    for lib in table[curdir]:
        own_name = "-l" + curdir.lower()
        if (src_test == 'test') or (lib != own_name):
            uniques[lib] = None
    for lib in uniques.keys():
        outstr += lib + '$(DIM)D' + ' '
    return outstr
