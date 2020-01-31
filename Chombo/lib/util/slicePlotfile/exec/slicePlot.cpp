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
#include "slice4Dplotfile.H"

using std::endl;

#include "UsingNamespace.H"

//using namespace Chombo;

void usageError(char* executableName)     
{ 
  std::cerr << " usage: " << executableName << " <filename> slice_dir1, slice_comp1 (slice_dir2 slice_comp2) \n"; exit(0); 
}


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

    Vector<int> slice_dir;
    Vector<int> slice_comp;
    if (argc < 4)
      {
        usageError(argv[0]);
      }

    char* in_file = argv[1];
    string fname4d(in_file);
    string slicefname(in_file);
    if (argc == 4)
      {
        slice_dir.resize(1);
        slice_comp.resize(1);
        slice_dir[0] = atoi(argv[2]);
        slice_comp[0] = atoi(argv[3]);
      }
    else if (argc == 6)
      {
        slice_dir.resize(2);
        slice_comp.resize(2);
        slice_dir[0] = atoi(argv[2]);
        slice_comp[0] = atoi(argv[3]);
        slice_dir[1] = atoi(argv[4]);
        slice_comp[1] = atoi(argv[5]);
      }
    else 
      {
        usageError(argv[0]);
      }
    

    // strip off "4d.hdf5" suffix and replace with either '3d.hdf5' or
    // "2d.hdf5"
    pout() << "4d plotfile in: " << fname4d.c_str() << endl;

    int position = slicefname.find("4d.hdf5");
    if (slice_dir.size() == 1)
      {
        slicefname.replace(position, 7, "3d.hdf5");
        pout() << "3d plotfile out: " << slicefname.c_str() << endl;
      } 
    else if (slice_dir.size() == 2)
      {
        slicefname.replace(position, 7, "2d.hdf5");
        pout() << "2d plotfile out: " << slicefname.c_str() << endl;
      }

    // call multidim function which does all the work...
    status = slice4Dplotfile(fname4d,
                             slicefname,
                             slice_dir, slice_comp);

  }  // end nested scope


#ifdef CH_MPI
  MPI_Finalize();
#endif

  return status;
}
