"""
This script attempts to generate a complete autoconf/automake/libtool-based
build system for any project that wants to use the Chombo library.  To use it,
type
  python make_example_metamakefiles.py <rootdir>
where rootdir is the root directory of your project.

Your project can contain C++ sources (with suffixes .cpp and .H) and
Chombo-Fortran sources (suffix .ChF and .fh).  The simplest sort of
project is a single directory with a single .cpp file.  But this
script can handle more complicated projects as well -- projects
distributed over several directories, like many of those under
Chombo's example directory.  Moreover, the build system we generate is
good for "make dist" as well as for out-of-source builds.

What we can't handle, though, are cases of:

  1. Multiple use of the same filename, e.g. dirA/Foo.H and dirB/Foo.H.
  2. Circular dependencies among directories, e.g. files in directory A
     #include'ing files in directory B and vice versa.
  3. Multiple definitions of a function, class, etc.
  4. Pretty much anything weirder than 1-3 above.
  5. Many cases where this script guesses wrong about the order in which
     directories need to be built -- see below.

If a directory contains a file or files with a "main(int argc...)" file,
we try to build an executable for each such file, and all other .cpp and .ChF
files in the same directory are compiled, and their .o files linked into
each executable.

If a directory does not contain any files with "main(int argc...)", we
roll all the .cpp and .ChF files in that directory into a library.  All
such libraries are linked into the executables mentioned in the paragraph
just above.

Coming now to #5 above: if we guess wrong about the order in which
directories need to be built, you will need to go into some of the
Makefile.am files and change the line that begins with "SUBDIRS".  And
then run the bootstrap script, followed by configure.

You will know the order is wrong if compilation fails because one of
your header files could not be found.  This will happen if a file in
directory A #include's a header from directory B -- counting on an
appropriate -I flag rather than specifying the relative path in the
#include directive -- and SUBDIRS lists directory A ahead of directory
B.

Another way SUBDIRS order will cause a build to fail is if the library
directories (the ones without main() functions) are not all built
before the directories that do have main() functions.  Once again, you
can fix this by changing SUBDIRS in the appropriate Makefile.am, and then
rerun bootstrap and configure.

A more intelligent program than this one would minimize such
subdirectory-order problems; someday perhaps this program will be like
that.  In the meantime, there are a few rules of thumb you can follow,
in the naming of your project's directories, that will increase the
chance of this script giving you a working build system:

  1. Have just one library directory (i.e. a directory with no main()
     functions), call it "src", and put it right under the root of
     your tree.
  2. If you must have more than one library directory, call the others
     src<something> (e.g. srcBob, srcEdna) and make sure they don't
     #include one another's headers (but they can #include headers from
     the src directory).
"""

import sys
import os
import re
import glob
import string


class WalkExtra:
    def __init__(self, indent='', all_lib_names=[] ):
        self.lib_dirs = []
        self.main_dirs = []
        self.indent = indent
        self.all_lib_names = all_lib_names


def walk1Deep( dirname, visit_func, extra_stuff ):
    apply( visit_func, ( extra_stuff, dirname, os.listdir(dirname) ) )

def hasMainFunc( filename ):
    """
    Returns True if file contains something that looks like a main() or
    main( int argc, char argv ) function.
    """
    main_regex = re.compile('.*main(.*argc.*argv).*')
    f = open(filename)
    for line in f.readlines():
        if main_regex.search( line ):
            return True
    return False
    

def makeMakefileAm( dirname, subdirs, main_cpps, lib_cpps, nonmain_nonlib_cpps,
                    lib_ChFs, depth ):
    makefile_am = open( dirname + '/Makefile.am', 'w' )

    # Directories names "src" or "src<something>" are good bets for
    # example-library code and should therefore be built first.
    if len(subdirs) > 0:
        def mysort(a,b):
            if a == 'src':
              return -1
            elif b == 'src':
              return 1
            elif a[:3] == 'src':
              return -1
            elif b[:3] == 'src':
              return 1
            return cmp(a,b)
        ordered_subdirs = subdirs[:]
        uniqer = {} # To remove duplicates...
        for s in ordered_subdirs : uniqer[s] = 0
        ordered_subdirs = uniqer.keys()
        ordered_subdirs.sort(mysort)
        makefile_am.write( 'SUBDIRS = ' + string.join( ordered_subdirs ) + '\n' )

    if len(lib_ChFs) + len(lib_cpps) + len(main_cpps) > 0:
        makefile_am.write( 'include ' + '../'*depth + 'Automake.rules\n' )

    if len(lib_ChFs) > 0:
        makefile_am.write( 'nodist_fort_HEADERS = ' )
        for f in lib_ChFs:
            bfname = f[:f.index('.ChF')]
            makefile_am.write( '\\\n    ' + bfname + '_F.H' )
        makefile_am.write( '\n' )
        makefile_am.write( 'fortdir = $(pkgincludedir)\n' )

        makefile_am.write( 'GENERATED_FORTRAN = ' )
        for f in lib_ChFs:
            bfname = f[:f.index('.ChF')]
            makefile_am.write( '\\\n    ' + bfname + '.f' )
        makefile_am.write( '\n' )
        makefile_am.write( 'EXTRA_DIST += *.ChF\n' )

    if len(main_cpps) > 0:
        makefile_am.write( 'bin_PROGRAMS = ' )
        for f in main_cpps:
            bfname = f[:f.index('.cpp')]
            makefile_am.write( ' ' + bfname )
        makefile_am.write('\n')

        for f in main_cpps:
            bfname = f[:f.index('.cpp')]
            makefile_am.write( bfname + '_SOURCES = ' )
            if len(lib_ChFs) > 0:
                makefile_am.write( '$(GENERATED_FORTRAN) ' )
            makefile_am.write( bfname + '.cpp' )
            for g in nonmain_nonlib_cpps:
                bgname = g[:g.index('.cpp')]
                makefile_am.write( ' ' + bgname + '.cpp' )
            makefile_am.write( '\n' )
        makefile_am.write( 'AM_LDFLAGS += -L$(pkglibdir) \n' )


    if len(lib_cpps) > 0:
        normalized_dirname = makeLTLibName( dirname )
        makefile_am.write( normalized_dirname + '_LTLIBRARIES = '
                         + 'lib' + normalized_dirname + '.la\n' )
        makefile_am.write( 'nodist_lib' + normalized_dirname + '_la_SOURCES = ' )
        if len(lib_ChFs) > 0:
            makefile_am.write( '$(GENERATED_FORTRAN)' )
        for f in lib_cpps:
            makefile_am.write( '\\\n    ' + f )
        makefile_am.write( '\n' )
        makefile_am.write( normalized_dirname + 'dir = $(pkglibdir)\n' )
        makefile_am.write( 'EXTRA_DIST += *.cpp\n' )

    if (len(lib_cpps) > 0) and (len(main_cpps) > 0):
        sys.stderr.write( "Warning: directory " + dirname + " contains both "
            + "cpp files with main() and cpp files without main()\n" )


    has_H = False
    if glob.glob( dirname + '/*.H' ) != []:
        has_H = True
        makefile_am.write( 'headers_HEADERS = $(srcdir)/*.H\n' )
        makefile_am.write( 'headersdir = $(pkgincludedir)\n' )
    if glob.glob( dirname + '/*.fh' ) != []:
        if has_H:
            makefile_am.write( 'headers_HEADERS += $(srcdir)/*.fh\n' )
        else:
            makefile_am.write( 'headers_HEADERS = $(srcdir)/*.fh\n' )
            makefile_am.write( 'headersdir = $(pkgincludedir)\n' )


def makeLTLibName( dirname ):
    return string.lstrip( dirname.replace( '/', '_' ), '._' )


def describeDirStructure( extra, dirname, files ):
    """
    Walk the tree under arg rootdir and classify the directories as to
    whether they are source (used for building a library) or exec (contain
    a cpp file with main() function in it).  A directory that has both is
    not considered a library directory; the source files that don't have
    a main() are simply considered binSOURCE's for the files that do.
    """
    bname = os.path.basename( dirname )
    if bname == 'CVS':
        return

    subdirs = filter( lambda f: os.path.isdir(dirname+'/'+f), files )

    has_cpp = bool(filter( lambda f: f.find('.cpp') != -1, files ))
    has_ChF = bool(filter( lambda f: f.find('.ChF') != -1, files ))

    has_lib_src=False
    has_lib_cpp=False
    has_main_cpp=False
    main_cpps = []
    nonmain_nonlib_cpps = []  # cpp's in dir that has a main() cpp.
    lib_cpps = []
    lib_ChFs = []
    if has_cpp:
        for cpp in glob.glob( dirname + '/*.cpp' ): 
            if hasMainFunc( cpp ):
                has_main_cpp = True
                main_cpps.append( os.path.basename(cpp) )
        for cpp in glob.glob( dirname + '/*.cpp' ): 
            if not hasMainFunc( cpp ):
                if has_main_cpp:
                    nonmain_nonlib_cpps.append( os.path.basename(cpp) )
                else:
                    has_lib_cpp = True
                    lib_cpps.append( os.path.basename(cpp) )
        if has_lib_cpp:
            extra.all_lib_names.append( makeLTLibName( dirname ) )

    if has_ChF:
        for chf in glob.glob( dirname + '/*.ChF' ):
            lib_ChFs.append( os.path.basename( chf ) )

    if has_main_cpp:
        extra.main_dirs.append( bname )
    if has_lib_cpp or has_ChF:
        has_lib_src = True
        extra.lib_dirs.append( bname )

    subdir_extra = WalkExtra( indent = extra.indent + '  ',
                              all_lib_names = extra.all_lib_names )
    for sd in subdirs:
        if sd == 'CVS': continue
        #print extra.indent, "Entering subdir", dirname+'/'+sd, "..."
        walk1Deep( dirname+'/'+sd, describeDirStructure, subdir_extra )
    if bool(subdir_extra.lib_dirs) or bool(subdir_extra.main_dirs):
        #print extra.indent, "lib subdirs of", bname, ":", subdir_extra.lib_dirs
        #print extra.indent, "main subdirs of", bname, ":", subdir_extra.main_dirs
        if not bname in extra.lib_dirs: # Source code >1 level down
            extra.lib_dirs.append( bname )

    if( len(main_cpps) + len(lib_cpps) + len(nonmain_nonlib_cpps) + len(lib_ChFs)
    +   len(subdir_extra.lib_dirs) + len(subdir_extra.main_dirs) > 0):
        makeMakefileAm( dirname, subdir_extra.lib_dirs + subdir_extra.main_dirs,
                        main_cpps, lib_cpps, nonmain_nonlib_cpps, lib_ChFs,
                        depth=len(extra.indent)/2 )


def fixupMakefileAms( all_lib_names, dirname, files ):
    """
    Every time you find a Makefile.am, look for "bin_PROGRAMS".  If it's there,
    then for every listed program, make an LDADD line that lists arg
    all_lib_names.
    """
    if 'Makefile.am' in files:
        ro = open( dirname + '/Makefile.am' )
        def grep( str, line ):
            return line.find(str) != -1
        bin_progs_line = filter( lambda line: grep('bin_PROGRAMS',line), ro.readlines() )
        if len(bin_progs_line) > 0:
            m = open( dirname + '/Makefile.am', 'a' )
            for bin in bin_progs_line[0].split()[2:]:
                m.write( bin + '_LDADD = -L$(CHOMBO_INSTALLDIR)/lib/Chombo $(LIBSRC_LIBS) ' )
                for lib in all_lib_names:
                    m.write( '-l' + lib + ' ' )
                    m.write( '-lg2c' )
                m.write( '\n' )
        

def findAllMakefileAms( topdir ):
    def visit_func( extra, dirname, files ):
        if 'Makefile.am' in files:
            extra['dirs_with_makefileam'
                 ].append( string.lstrip(
                    dirname[len(extra['topdir'])+1:] + '/Makefile', '/'))
    extra = {'dirs_with_makefileam':[], 'topdir':topdir}
    os.path.walk( topdir, visit_func, extra )
    return extra['dirs_with_makefileam']


def  makeConfigureIn( topdir ):
    """
    Generate the configure.in, starting from a template and just adding
    AC_CONFIG_FILES and AC_OUTPUT lines.

    It's assumed this script is in the same directory as configure.pre; they should
    both be at the root of the Chombo source tree, and in share/Chombo in any
    Chombo install tree.
    """
    confpre = os.path.dirname(sys.argv[0]) + '/../../../configure.pre'
    if confpre[0] == '/':
        confpre = confpre[1:]
    infile  = open( confpre )
    outfile = open( topdir + '/configure.in', 'w' )

    in_lines = infile.readlines()
    curline = 0
    while True:  # No do...while in Python.
        outfile.write( in_lines[curline] )
        curline += 1
        if in_lines[curline-1][0:7] == "AC_INIT":
            break

    #
    # Users of the configure.in we're writing will need to tell it where their
    # Chombo install tree is.  That's where Make.defs and transformation scripts are.
    #
    outfile.write( 'if test x$CHOMBO_INSTALLDIR = x ; then\n' )
    outfile.write( '  echo " no CHOMBO_INSTALLDIR" \n' )
    outfile.write( 'fi\n' )

    outfile.write( 'if test ! -f $CHOMBO_INSTALLDIR/share/Chombo/Make.defs ; then\n' )
    outfile.write( '  echo " no $CHOMBO_INSTALLDIR/share/Chombo/Make.defs" \n' )
    outfile.write( 'fi\n' )


    outfile.write( 'if test x$CHOMBO_INSTALLDIR = x -o ! -f $CHOMBO_INSTALLDIR/share/Chombo/Make.defs ; then\n' )
    outfile.write( '  echo "*****************************************" \n' )
    outfile.write( '  echo "Error: you must pass configure a definition of CHOMBO_INSTALLDIR"\n' )
    outfile.write( '  echo "and it must indicate the root of your Chombo install tree" \n' )
    outfile.write( '  echo "e.g. \'./configure CHOMBO_INSTALLDIR=\$HOME/Chombo/install\'"\n' )
    outfile.write( '  echo ""\n'                                                            )
    outfile.write( '  echo "If you think this message is in error, check that under your"\n')
    outfile.write( '  echo "CHOMBO_INSTALLDIR you have a file called Make.defs.  If you"\n' )
    outfile.write( '  echo "do not (but, say, you do seem to have some of the libraries"\n')
    outfile.write( '  echo "and header files), then it is possible your Chombo build just"\n' )
    outfile.write( '  echo "did not run to completion."\n' )
    outfile.write( '  echo "*****************************************" \n' )
    outfile.write( '  exit 1 \n' )
    outfile.write( 'fi\n' )


    #
    # Substitute something for the project name.
    #
    while in_lines[curline][0:16] != 'AM_INIT_AUTOMAKE':
        outfile.write( in_lines[curline] )
        curline += 1
    outfile.write( 'AM_INIT_AUTOMAKE('+os.path.basename(topdir)+', 0.1.0 )\n' )
    curline += 1

    #
    # Throw away configure.pre lines that control what is and isn't in a
    # "small build".
    #
    for line in in_lines[curline:]:
        if line[0:11] == "#SMALLBUILD"  or  line[0:12] == "#!SMALLBUILD" :
            continue
        outfile.write( line )

    #
    # Write out the paths to the Makefiles we want generated.
    #
    outfile.write( 'AC_CONFIG_FILES(\n' )
    makefile_ams = findAllMakefileAms( topdir )
    for m in makefile_ams:
        outfile.write( '    ' + m + '\n' )
    outfile.write( ')\n' )
    outfile.write( 'AC_OUTPUT\n' )
    
            
if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.stderr.write( "**********************************************************\n")
        sys.stderr.write( "Usage: python make_example_metamakefiles.py example_dir   \n")
        sys.stderr.write( "  where example_dir is the name of the directory in which \n")
        sys.stderr.write( "  you want to generate an autoconf build system.\n")
        sys.stderr.write( "\n" )
        sys.stderr.write( "It's important that when you run this script, your working\n")
        sys.stderr.write( "directory be one above example_dir.\n")
        sys.stderr.write( "**********************************************************\n")
        sys.exit(1)

    topdir = string.rstrip( sys.argv[1], '/' )
    chombodir = os.path.dirname(sys.argv[0]) + '/../../..'
    print "chombodir=", chombodir
    walk_extra = WalkExtra( indent='', all_lib_names=[] )
    walk1Deep( topdir, describeDirStructure, walk_extra )

    # Now that you know the full set of LT libraries, fix up the Makefile.am's
    # of the executables, so they all link to all those LT libraries.
    os.path.walk( topdir, fixupMakefileAms, walk_extra.all_lib_names )

    # Generate the configure.in, starting from a template and just adding
    # AC_CONFIG_FILES and AC_OUTPUT lines.
    makeConfigureIn( topdir )

    assert( os.path.exists( chombodir + '/lib/mk/autoconf/Automake.rules' ) )
    os.system( "sed 's/makefile\.anag//' " + chombodir + "/lib/mk/autoconf/Automake.rules > "
                + topdir + "/Automake.rules" )
    os.system( "cp " + chombodir + "/lib/mk/autoconf/bootstrap-generated-example.sh "
             + topdir + "/bootstrap\n" )
    os.system( "cp " + chombodir + "/lib/mk/autoconf/zap-generated-example.sh "
             + topdir + "/zap\n" )

    # Now that you've created configure.in and all the Makefile.am's, run
    # GNU autotools on them.
#   os.system( 'cd ' +  topdir + '; ./bootstrap' )
