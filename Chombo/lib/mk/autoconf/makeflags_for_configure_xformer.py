"""
Clean up $(MAKEFLAGS) to make them palatable as options to the configure
script.  In particular, this means removing anything that's not a key=value
pair.
"""

import sys
import re
import string

regex = re.compile( '[A-Za-z_][A-Za-z0-9_]* *= *[A-Za-z0-9_][A-Za-z0-9_]*' )

line = sys.stdin.readline()

refindall = re.findall( regex, line )
key_value_pairs = string.join( refindall )
print "key_value_pairs=", key_value_pairs

# Now separate out the rest of the stuff
remainder = line[:]
all_else = ''
for tok in refindall:
    all_else += remainder[:remainder.index(tok)]
    remainder = remainder[remainder.index(tok) + len(tok):]
all_else += remainder
print "all_else(final):", all_else
