#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>
#include <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "SPMD.H"
#include "UGIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "Misc.H"
#include "Vector.H"
#include "REAL.H"
#include "Box.H"
#include "BoxIterator.H"

#include "UsingNamespace.H"

Real getDataVal(const IntVect& a_iv)
{
  Real retval  = 7.23;
  Real dx = 0.001;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Real arg = Real(idir + a_iv[idir]*a_iv[idir]);
      retval += sin(dx*arg)+ 2.*cos(dx*arg);
    }
  return retval;
}
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  {//scoping trick

    //set grid spacing,
    // number of points in each direction, domain length
    int nx = 64;
    // number of processors
    int nproc = numProc();
    // make maximum box size
    int maxsize = Max(nx/(2*nproc), 4);
    // this is the domain of the computation
    Box domain(IntVect::Zero, (nx-1)*IntVect::Unit);
    // this is the list of boxes in the layout
    Vector<Box> vbox;
    int blockfac = 1;
    domainSplit(domain, vbox, maxsize, blockfac);
    //load balance the boxes
    Vector<int> procAssign;
    LoadBalance(procAssign, vbox);
    /// create layout
    DisjointBoxLayout dbl(vbox, procAssign);
    ///make the data to output and set its values.
    LevelData<FArrayBox> data(dbl, 1);
    DataIterator dit = dbl.dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
      {
        FArrayBox& fab = data[dit()];
        const Box& fabbox = fab.box();
        BoxIterator bit(fabbox);
        for (bit.reset(); bit.ok(); ++bit)
          {
            const IntVect& iv = bit();
            fab(iv, 0) = getDataVal(iv);
          }
      }
#ifdef CH_USE_HDF5
    //output data to a file
    string filename("dataout.hdf5");
    WriteUGHDF5(filename, dbl, data, domain);

    //Read it back in.
    DisjointBoxLayout dblin;
    LevelData<FArrayBox> datain;
    Box domainin;
    ReadUGHDF5(filename, dblin, datain, domainin);
    //check to see that it matches.
    if (domainin != domain)
      {
        cerr << "domains do not match" << endl;
        return -1;
      }
    if (datain.nComp() != 1)
      {
        cerr << "input data has the wrong number of components" << endl;
        return -2;

      }
    DataIterator ditin = dblin.dataIterator();
    for (ditin.reset(); ditin.ok(); ++ditin)
      {
        const FArrayBox& fabin = datain[ditin()];
        const Box& fabbox = fabin.box();
        BoxIterator bit(fabbox);
        for (bit.reset(); bit.ok(); ++bit)
          {
            const IntVect& iv = bit();
            Real rightans = getDataVal(iv);
            Real dataans = fabin(iv, 0);
            Real eps = 1.0e-9;
            if (Abs(dataans - rightans) > eps)
              {
                cerr << "data does not match" << endl;
                return -3;
              }
          }
      }
#endif // CH_USE_HDF5
  }//end scoping trick
  cout << "single level IO seems to work" << endl;
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return(0);
}
