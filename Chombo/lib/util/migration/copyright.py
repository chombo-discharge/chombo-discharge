#!/usr/bin/env python

#      _______              __
#     / ___/ /  ___  __ _  / /  ___
#    / /__/ _ \/ _ \/  V \/ _ \/ _ \
#    \___/_//_/\___/_/_/_/_.__/\___/
#    Please refer to Copyright.txt, in Chombo's root directory.


#
# Strips out all our diverse copyright notices and replaces them with a
# single, standardized one -- the one you see five lines up from here in
# fact.
#
# Works on all Chombo .cpp, .H and .ChF files, except a couple under example
# where it chomps off a few #include lines (and I fixed those by hand).
#
# The "main" of this program takes a single argument -- the name of a file.
# The file is not modified, but a new one, with the proper new-style copyright
# notice is created: its name is the same as the original, plus ".new".
# If you want to change all the files in a tree, a sh script like this will
# work (after you put this file -- copyright.py -- on your PATH):
"""
#!/bin/sh
for ext in H cpp ChF ; do
  for f in `find . -name "*.$ext" | grep -v ccse` ; do
    copyright.py $f
    mv ${f}.new $f
  done
done
"""
# You'll get a few warnings along the way.  The ones about files that didn't
# have any copyright notice to begin with are harmless (and they're the only
# ones I saw).

import sys
import string

g_boilerplate = []
g_boilerplate.append(
"""\
// CHOMBO Copyright (c) 2000-2004, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory (subject to
// receipt of any required approvals from U.S. Dept. of Energy).  All
// rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// (1) Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// (2) Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
// (3) Neither the name of Lawrence Berkeley National Laboratory, U.S.
// Dept. of Energy nor the names of its contributors may be used to endorse
// or promote products derived from this software without specific prior
// written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// You are under no obligation whatsoever to provide any bug fixes,
// patches, or upgrades to the features, functionality or performance of
// the source code ("Enhancements") to anyone; however, if you choose to
// make your Enhancements available either publicly, or directly to
// Lawrence Berkeley National Laboratory, without imposing a separate
// written license agreement for such Enhancements, then you hereby grant
// the following license: a non-exclusive, royalty-free perpetual license
// to install, use, modify, prepare derivative works, incorporate into
// other computer software, distribute, and sublicense such Enhancements or
// derivative works thereof, in binary and source code form.
//
// TRADEMARKS. Product and company names mentioned herein may be the
// trademarks of their respective owners.  Any rights not expressly granted
// herein are reserved.
//\
""".split('\n'))

g_boilerplate.append( g_boilerplate[-1][:] )
g_boilerplate[-1][0] = string.replace( g_boilerplate[-1][0], '2004', '2006' )

g_boilerplate.append(
"""\
// This software is copyright (C) by the Lawrence Berkeley
// National Laboratory.  Permission is granted to reproduce
// this software for non-commercial purposes provided that
// this notice is left intact.
// 
// It is acknowledged that the U.S. Government has rights to
// this software under Contract DE-AC03-765F00098 between
// the U.S.  Department of Energy and the University of
// California.
//
// This software is provided as a professional and academic
// contribution for joint exchange. Thus it is experimental,
// is provided ``as is'', with no warranties of any kind
// whatsoever, no support, no promise of updates, or printed
// documentation. By using this software, you acknowledge
// that the Lawrence Berkeley National Laboratory and
// Regents of the University of California shall have no
// liability with respect to the infringement of other
// copyrights by any part of this software.
//\
""".split('\n'))

g_boilerplate.append(
"""\
C  CHOMBO Copyright (c) 2000-2004, The Regents of the University of
C  California, through Lawrence Berkeley National Laboratory (subject to
C  receipt of any required approvals from U.S. Dept. of Energy).  All
C  rights reserved.
C
C  Redistribution and use in source and binary forms, with or without
C  modification, are permitted provided that the following conditions are met:
C
C  (1) Redistributions of source code must retain the above copyright
C  notice, this list of conditions and the following disclaimer.
C  (2) Redistributions in binary form must reproduce the above copyright
C  notice, this list of conditions and the following disclaimer in the
C  documentation and/or other materials provided with the distribution.
C  (3) Neither the name of Lawrence Berkeley National Laboratory, U.S.
C  Dept. of Energy nor the names of its contributors may be used to endorse
C  or promote products derived from this software without specific prior
C  written permission.
C
C  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
C  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
C  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
C  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
C  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
C  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
C  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
C  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
C  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
C  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C
C  You are under no obligation whatsoever to provide any bug fixes,
C  patches, or upgrades to the features, functionality or performance of
C  the source code ("Enhancements") to anyone; however, if you choose to
C  make your Enhancements available either publicly, or directly to
C  Lawrence Berkeley National Laboratory, without imposing a separate
C  written license agreement for such Enhancements, then you hereby grant
C  the following license: a non-exclusive, royalty-free perpetual license
C  to install, use, modify, prepare derivative works, incorporate into
C  other computer software, distribute, and sublicense such Enhancements or
C  derivative works thereof, in binary and source code form.
C
C  TRADEMARKS. Product and company names mentioned herein may be the
C  trademarks of their respective owners.  Any rights not expressly granted
C  herein are reserved.
C\
""".split('\n'))


g_boilerplate.append(
"""\
C     This software is copyright (C) by the Lawrence Berkeley
C     National Laboratory.  Permission is granted to reproduce
C     this software for non-commercial purposes provided that
C     this notice is left intact.
C
C     It is acknowledged that the U.S. Government has rights to
C     this software under Contract DE-AC03-765F00098 between
C     the U.S.  Department of Energy and the University of
C     California.
C
C     This software is provided as a professional and academic
C     contribution for joint exchange. Thus it is experimental,
C     is provided ``as is'', with no warranties of any kind
C     whatsoever, no support, no promise of updates, or printed
C     documentation. By using this software, you acknowledge
C     that the Lawrence Berkeley National Laboratory and
C     Regents of the University of California shall have no
C     liability with respect to the infringement of other
C     copyrights by any part of this software.
C\
""".split('\n'))


def logo( filename ):
    chombo = [
'#ifdef CH_LANG_CC',
'/*',
'*      _______              __',
'*     / ___/ /  ___  __ _  / /  ___',
'*    / /__/ _ \\/ _ \\/  V \\/ _ \\/ _ \\',
'*    \\___/_//_/\\___/_/_/_/_.__/\\___/',
"*    Please refer to Copyright.txt, in Chombo's root directory.",
'*/',
'#endif']

    result = ''
    if filename[-4:] == '.ChF':
        for i in (8,7,1,0):
            del chombo[i]

    for line in chombo:
        if filename[-4:] == '.ChF':
            line = string.replace(line,'*','C')
        result += line + '\n'
    return result


"""
Close-enough match between line in source file, and line we're looking for.
"""
def goodmatch( str1, str2 ):
    gstr1 = str1[1:].strip() # Ignores ^C vs ^!, and number blank spaces after C
    gstr2 = str2[1:].strip() # Ignores ^C vs ^!, and number blank spaces after C
    return gstr1 == gstr2


def stripBoilerplate( filename ):
    result = []
    boilerplate_line=0
    overall_line=-1
    bp_model = None # element of g_boilerplate -- the one found

    
    #
    # Insert pointer to copyright notice.
    #
    result.append( logo(filename) )


    #
    # Find which line copyright notice ends at, and copy file from that point
    # on.  If you can't find anything like a copyright notice, then copy the
    # entire file.
    #

    last_line_of_copyright = 0
    include_guards = []

    text = open(filename).readlines()
    for line in text:
        overall_line += 1

        # Don't lose the include guards, if they're at very top of file:
        if( (overall_line < 3) 
        and ( (-1 != line.find('#ifndef'))
           or (-1 != line.find('#define')))):
            include_guards.append( line )

        if not bp_model:
            for bp in g_boilerplate:
                if goodmatch( line[:-1], bp[0] ):
                    bp_model = bp
                    last_line_of_copyright = overall_line + len(bp_model)
                    break

            if overall_line == 20:
                print "Warning:(", filename, ") haven't found boilerplate yet."
                break

    # Go through once again, this time copying everything from the end of the
    # copyright notice.
    overall_line = -1
    if (last_line_of_copyright != 0) and (len(include_guards)==2):
        result += include_guards
    if -1 != text[last_line_of_copyright].find( '#endif' ): # matches CH_LANG_CC
        text[last_line_of_copyright] = '\n'                 # we removed.
    for iline in range(last_line_of_copyright, len(text)):
        result += text[iline]
        
    outfile = open( filename + '.new', 'w' )    
    for rl in result:
        outfile.write(rl)


if __name__ == '__main__':
    stripBoilerplate( sys.argv[1] )
