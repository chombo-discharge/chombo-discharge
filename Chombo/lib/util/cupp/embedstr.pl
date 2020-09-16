#!/usr/bin/perl

####  _______              __
#### / ___/ /  ___  __ _  / /  ___
####/ /__/ _ \/ _ \/  ' \/ _ \/ _ \
####\___/_//_/\___/_/_/_/_.__/\___/
####  Please refer to Copyright.txt, in Chombo's root directory.

# ----------------------------------------------------------------------
# Perl script to generate an embedded array of bytes from an input
# string and store them into a variable in a header file
#
# Usage:  perl embedstr.pl <inputfile> output_prefix
#   Note, the output_prefix can include a path but should not have a
#   '.' or suffix
#
# Use to embed Cuda executable byte strings.  See also fatbinary
# excutable distributed with Cuda and used with -no-asm.
# ----------------------------------------------------------------------

use strict;
use warnings;

my $wordBytes = 8;

if (@ARGV < 2)
{
  printf die "Usage: embedstr.pl <inputfile> output_prefix\n";
}

my $oname = $ARGV[1] . "_CUX.H";
open(my $ifh, "<:raw", $ARGV[0]) or die "embedstr.pl: Unable to open file $ARGV[0]\n";
open(my $ofh, ">$oname") or die "embedstr.pl: Unable to open file $oname\n";

my $tmpstr = $ARGV[1];
$tmpstr =~ m{([^/]+)$};  # Anchor at end, all characters not including /
my $varname = $1;
my $ucname = uc($varname) . "_CUX";
print {$ofh} "#ifndef _${ucname}_H_\n";
print {$ofh} "#define _${ucname}_H_\n";
print {$ofh} "\n// This is a binary module for the GPU\n";
print {$ofh} "static const unsigned long long $varname\[\] = {\n";
# First word
{
    read $ifh, (my $word), $wordBytes;
    my $uintword = unpack('Q', $word);
    printf {$ofh} "0x%016xull", $uintword;
}
# Remainder
my $wmod = 0;
while (read $ifh, (my $word), $wordBytes)
{
    my $len = length($word);
    if ($len <= 0)
    {
        die "embedstr.pl: Invalid word size\n";
    }
    while ($len < 8)
    {
        $word .= "\0";
        ++$len;
    }
    my $uintword = unpack('Q', $word);
    printf {$ofh} ",";
    if (++$wmod == 4)
    {
        printf {$ofh} "\n";
        $wmod = 0;
    }
    printf {$ofh} "0x%016xull", $uintword;
}
close($ifh);
# This would always ensure num_word % 4 == 0 ... but there is no good reason for
# doing it.
#if ($wmod != 3)
#{
#    while (++$wmod < 4)
#    {
#        printf {$ofh} ",%016xull", 0;
#    }
#}
print {$ofh} "\n};\n";
print {$ofh} "\n#endif\n";
close($ofh);
