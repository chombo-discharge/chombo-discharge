#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//===========================================================================
// injectPlot.cpp
//
//===========================================================================
#include <iostream>

#define CH_SPACEDIM 2


// includes
#define CH_SPACEDIM 2
#include "SPMD.H"
#include "AMRIO.H"
#include "CH_Attach.H"
#include "convert1DplotTo2D.H"

using std::endl;

#include "UsingNamespace.H"

//using namespace Chombo;


//===========================================================================
// simple code to read in 1D plotfile, inject to 2D, then write 
// 2d plotfile
//
//===========================================================================
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  int status = 0;

  //registerDebugger();

  { // Begin nested scope

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
    int rank, number_procs;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
    rank=0;
    number_procs=1;
#endif

    int num_width=1;
    if (argc < 2)
      { std::cerr << " usage: " << argv[0] << " <filename>, (<num_width>)\n"; exit(0); }
    char* in_file = argv[1];
    string fname1d(in_file);
    string fname2d(in_file);

    if (argc > 2)
      {
        num_width = atoi(argv[2]);
      }
    if (num_width != 1)
      {
        MayDay::Error("num_width != 1 not implemented yet");
      }

    // strip off "1d.hdf5" suffix and replace with '2d.hdf5'
    int position = fname2d.find("1d.hdf5");
    fname2d.replace(position, 7, "2d.hdf5");

    pout() << "1d plotfile in: " << fname1d.c_str() << endl;
    pout() << "2d plotfile out: " << fname2d.c_str() << endl;
    pout() << "2d width = " << num_width << endl;

    // call multidim function which does all the work...
    status = convert1DPlotTo2D(fname1d,
                                   fname2d,
                                   num_width);

  }  // end nested scope


#ifdef CH_MPI
  MPI_Finalize();
#endif

  return status;
}
