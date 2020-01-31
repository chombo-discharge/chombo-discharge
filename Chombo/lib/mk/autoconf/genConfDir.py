import sys
import string

MAKEFLAGS = sys.argv[1:]
MAKEFLAGS.sort() # To get things in a canonical order
result = ''
for tok in MAKEFLAGS:
    result += string.replace(tok, '=', '_') + '.'
if MAKEFLAGS == []:
    result = 'defaults.'
print result[:-1]
