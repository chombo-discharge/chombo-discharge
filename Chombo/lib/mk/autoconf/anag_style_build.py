import re
import sys
import os
import string


def doingMultidim():
    """
    Returns True if MULTIDIM is set to TRUE either in Make.defs or in $(MAKEFLAGS).
    """
    if( ('DIM' in makeflags()['varvaluedict'].keys())
    and (makeflags()['varvaluedict']['DIM']=='TRUE') ):
        sys.stderr.write( "Error: If MULTIDIM==TRUE, DIM is forbidden " +
                          "on the make command line.\n" )
        sys.exit(1)

    return configOptions()['MULTIDIM'] == 'TRUE'


def thisConfigDir():
    """
    Returns the name of the (out-of-source) build directory in which to make this
    configuration.  That name is composed from the options, listed in Make.defs,
    that were overridden on the "gmake -f makefile.anag" command line.
    """

    defaults = defaultConfigOptions()
    sorted_keys = defaults.keys()
    sorted_keys.sort() # ensures a canonical order
    value_MAKEFLAGS = makeflags()['varvaluedict']

    result = ''
    for k in sorted_keys:
        if (k in value_MAKEFLAGS.keys()) and (value_MAKEFLAGS[k] != defaults[k]):
            result += '.' + k + '_' + value_MAKEFLAGS[k]
    if result == '':
        result = '.defaults'
    return result[1:]  # Strips off leading '.'


def defaultConfigOptions():
    """
    Return a dictionary showing all the configuration options as given in Make.defs.
    """
    defaults = {}
    junklineregexp = re.compile( '^ *\t*#.*$|^ *\t*$' )
    for default in open( topSrcDir() + '/lib/mk/autoconf/Make.defs' ).readlines():
        if junklineregexp.search( default ):
            continue
        temp = default.split('=')
        defaults[ temp[0] ] = temp[1][:-1]
    return defaults


def makeflags():
    """
    This script is launched by makefile.anag, which passes it the following sys.argv:
    '$(MAKE) ||| $(MAKEFLAGS) target'.  This function parses out various parts of that,
    returning a dictionary with these keys:
      'target': target
      'make': original $(MAKE)
      'makeflags': original $(MAKEFLAGS)
      # ...and these parts of makeflags...
      'varvaluedict': a dictionary showing all KEY=VALUE (e.g. DEBUG=FALSE)
                      configuration options given in $(MAKEFLAGS);
      'varvaluestr': the KEY=VALUE part of $(MAKEFLAGS), as a string.
      'other': everything else (from $(MAKEFLAGS), as a string, for example if make
               was started with -j4, 'other' will have something like
               "--jobserver-fds=3,5 -j".
    """
    result = {}

    result['target'] = sys.argv[-1:][0]
    assert( result['target'] in ('all', 'install', 'clean', 'uninstall', 'realclean') )

    separator_pos = sys.argv.index('|||')
    assert( separator_pos == 2 )

    result['make'] = sys.argv[1]

    MAKEFLAGS = string.join(sys.argv[separator_pos+1:-1])
    result['makeflags'] = MAKEFLAGS

    regex = re.compile( '[A-Z][A-Z0-9_]* *= *[A-Za-z0-9_]+' )
    refindall = re.findall( regex, MAKEFLAGS )
    result['varvaluestr'] = string.join(refindall)

    result['varvaluedict'] = {}
    if refindall:
        for arg in refindall:
            temp = arg.split('=')
            result['varvaluedict'][temp[0]] = temp[1]

    # Now separate out the rest of the stuff
    remainder = MAKEFLAGS[:]
    result['other'] = ''
    for tok in refindall:
        result['other'] += remainder[:remainder.index(tok)]
        remainder = remainder[remainder.index(tok) + len(tok):]
    result['other'] += remainder
    return result


def configOptions():
    """
    Return a dictionary showing all the configuration options.  The keys are all the
    items in Make.defs to the left of the '=' signs.  The values are the ones
    in Make.defs, except when overridden in $(MAKEFLAGS).
    """
    result = defaultConfigOptions()
    for k in makeflags()['varvaluedict'].keys():
        result[k] = makeflags()['varvaluedict'][k]
    return result


def topSrcDir():
    return os.path.abspath(os.path.dirname(sys.argv[0]) + '/../../..')
    # (We're in lib/mk/autoconf)


def thisSubDir():
    """
    Returns the name, relative to topsrcdir, of the directory from which this script
    was called.
    """
    top_src_dir = topSrcDir()
    cwdpth = os.getcwd()
    if top_src_dir == cwdpth:
        result = ''
    else:
        result = cwdpth[len(top_src_dir)+1:]
    return result


def multidimCodeBelow( dirname, multidim_dirs ):
    """
    Return True if in this directory or below it on the tree there's any multidim
    code, which is defined by the presence of directories named in the tuple argument
    multidim_dirs.
    """
    foundit = [False] # Has to be a structure, so it gets passed by reference.

    if os.path.basename(dirname) in multidim_dirs:
        return True

    def visit( arg, curdir, names ):
        if os.path.basename(curdir) in multidim_dirs:
            arg[0] = True
            return

    os.path.walk( dirname, visit, foundit )
    return foundit[0]

def unidimCodeBelow( dirname, multidim_dirs, ignore_dirs ):
    """
    Return True if in this directory or below it on the tree there's anything *but*
    multidim code, which is defined by the presence of directories named in the tuple
    argument multidim_dirs (but excluding directories names in the tuple argument
    ignore_dirs).
    """
    foundit = [False] # Has to be a structure, so it gets passed by reference.

    if( (not os.path.basename(dirname) in multidim_dirs)
    and (not os.path.basename(dirname) in ignore_dirs) ):
        return True

    def visit( arg, curdir, names ):
        if( (not os.path.basename(curdir) in multidim_dirs)
        and (not os.path.basename(curdir) in ignore_dirs) ):
            arg[0] = True
            return

    os.path.walk( dirname, visit, foundit )
    return foundit[0]


if __name__ == '__main__':
    """
    argv[1:] is $(MAKE) "|||" $(MAKEFLAGS) target.
    """
    buildsdir = topSrcDir() + '/builds'
    thisconfdir = buildsdir + '/' + thisConfigDir()

    if makeflags()['target'] == 'realclean':
        os.system( 'rm -rf ' + thisconfdir )
        sys.exit(0)

    # Create root builds dir, if necessary.
    if not os.path.isdir( buildsdir ):
        sys.stderr.write( "builds directory does not exist.  Making it now...\n" )
        os.mkdir( buildsdir )

    spacedim_range =  range( int(defaultConfigOptions()['CH_MIN_SPACEDIM']),
                             int(defaultConfigOptions()['CH_MAX_SPACEDIM'])+1 )
    
    # Create and configure builds subdir for this configuration, if
    # necessary.
    if not os.path.isdir( thisconfdir ):
        os.mkdir( thisconfdir )
        if doingMultidim():
            multidimdir = thisconfdir + '/multidim'
            os.system( 'mkdir -p ' + multidimdir )
            os.system( 'cd ' + multidimdir + ';' +
                   topSrcDir() + '/configure --prefix=`pwd`/../install' +' '
                      + makeflags()['varvaluestr'] + ' DIM=0')                
            for dim in spacedim_range:
                dimdir = thisconfdir + '/d' + str(dim)
                os.system( 'mkdir -p ' + dimdir )
                os.system( 'cd ' + dimdir + ';' +
                       topSrcDir() + '/configure --prefix=`pwd`/../install '
                          + '--without-tests --without-examples '
                          + makeflags()['varvaluestr'] + ' DIM=' + str(dim))
        else:
            os.system( 'cd ' + thisconfdir + ';' +
                       topSrcDir() + '/configure --prefix=`pwd`/install' + ' '
                          + makeflags()['varvaluestr'] )
    # Make.
    if doingMultidim():
        multidim_dirs = ('MultiDim', 'smallMultiDim') # Add to this, as necessary.
        ignore_dirs = ('CVS')
        if unidimCodeBelow(os.getcwd(), multidim_dirs, ignore_dirs):
            sys.stderr.write( "unidim code found below " + thisSubDir() + '\n' )
            for dim in spacedim_range:
                retval = \
                  os.system( 'cd '+ thisconfdir + '/d' + str(dim) + '/' + thisSubDir() + ';'
                    + makeflags()['make'] + ' ' + makeflags()['makeflags'] + ' '
                    + makeflags()['target'])
                if retval != 0:
                    break
        if multidimCodeBelow(os.getcwd(), multidim_dirs):
            sys.stderr.write( "multidim code found below " + thisSubDir() + '\n' )
            os.system( 'cd ' + thisconfdir + '/multidim/' + thisSubDir() + ';'
                + makeflags()['make'] + ' ' + makeflags()['makeflags'] + ' '
                + makeflags()['target'])
    else:    
        os.system( 'cd ' + thisconfdir + '/' + thisSubDir()
                  + '; ' + makeflags()['make'] + ' ' + makeflags()['makeflags'] + ' '
                  + makeflags()['target'] )
