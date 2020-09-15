import sys
import string

def stripQuotes( file_lines, quote_chars, mirror=None, replace=False ):
    """
    Strips out everything between pairs of quote_chars.  Arg quote_chars is
    typically ', ", or \"\"\".

    Arg file_lines is what you'd get from a file.readlines() on the file.

    If arg mirror==1, then the close-quote is taken as quote_chars.reverse()
    (useful for stripping out /* ... */ pairs).

    Returns the lines of the file -- minus the material between matching quotes
    -- unless optional arg replace is True, in which case everything but '\n'
    is converted to blanks.
    """
    result = ''  # Will be a list of strings -- lines of the file.

    # We'll work on file_lines joined into a single long string, as that
    # eliminates the complications of finding pairs of quotes that span
    # multiple lines.
    # After we're done, we'll break the result back up at the '\n's, so the
    # other filters will work properly.
    joined_lines = string.join( file_lines, '' )

    if mirror == 1:
        reversed_quote_chars = list(quote_chars)
        reversed_quote_chars.reverse()
        reversed_quote_chars = string.join( reversed_quote_chars, '' )

    open_pos = joined_lines.find( quote_chars )
    if open_pos == -1:
        # Didn't find any quote chars.  No problem...
        return file_lines
    result = result + joined_lines[:open_pos]

    while 1:
        if mirror == 1:
            close_pos = joined_lines.find( reversed_quote_chars,
                                           open_pos + len(quote_chars))
        else:
            close_pos = joined_lines.find( quote_chars,
                                           open_pos + len(quote_chars))
        if close_pos == -1:
            sys.stderr.write( "!!! Failed to find closing quote near character"
                + str(open_pos) + ". Context was"
                + joined_lines[open_pos:open_pos+100]
                + ". mirror=" + str(mirror) + " quote_chars=" + quote_chars)
            return result

        if replace == True:
            result += blankChars(joined_lines[open_pos:
                                              close_pos+len(quote_chars)])

        open_pos = joined_lines.find( quote_chars, close_pos + len(quote_chars))
        if open_pos == -1:
            result += joined_lines[close_pos+len(quote_chars):]
            break
        else:
            result += joined_lines[close_pos+len(quote_chars):open_pos]

    # Convert back to a list of '\n'-terminated strings.
    result = result.split('\n')
    for i in range(0, len(result)):
        result[i] = result[i] + '\n'
    return result


if __name__ == '__main__':
    """
    Reads from stdin, strips off /**/-style comments, and dumps to stdout.
    """
    clean = stripQuotes( sys.stdin.readlines(), '/*', mirror=1 );
    for line in clean:
        print line[:-1]
