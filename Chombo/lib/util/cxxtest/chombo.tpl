// -*- C++ -*-
// This custom Chombo test runner using CxxTest template files.
// It makes the command-line arguments available for MPI
// initialization and passing in of an input-file name.

#include <stdlib.h>
#include <cxxtest/ChomboPrinter.h>

int* _argcptr = NULL;
char*** _argvptr = NULL;

int main( int argc, char *argv[] )
{
  _argcptr = &argc;
  _argvptr = &argv;
  return CxxTest::ChomboPrinter().run();
}

// The CxxTest "world"
<CxxTest world>

