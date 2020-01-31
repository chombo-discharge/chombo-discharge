#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include  <iostream>
using std::cout;
using std::endl;

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
#include "heatfortF_F.H"

#include "UsingNamespace.H"

int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  {//scoping trick

    // number of points in each direction, domain length
    //and diffusion coefficient
    int nx = 64;
    Real domainLen = 1.0;
    Real coeff = 1.0e-3;
    //define stopping conditions
    Real tfinal = 3.33;
    int nstepmax = 100;
    //set grid spacing and time step
    Real dx = domainLen/nx;
    Real dt = 0.8*dx*dx/(2.*SpaceDim*coeff);
    // number of processors
    int nproc = numProc();
    // make maximum box size
    int maxsize = Max(nx/nproc, 4);
    // this is the domain of the computation
    Box domain(IntVect::Zero, (nx-1)*IntVect::Unit);
    // this is the list of boxes in the layout
    Vector<Box> vbox;
    domainSplit(domain, vbox, maxsize);
    //load balance the boxes
    Vector<int> procAssign;
    LoadBalance(procAssign, vbox);
    /// create layout
    DisjointBoxLayout dbl(vbox, procAssign);
    ///make the data with one ghost cell for convenience
    LevelData<FArrayBox> phi(dbl, 1, IntVect::Unit);
    LevelData<FArrayBox> lph(dbl, 1, IntVect::Zero);

    //set phi to 1 for an initial condition
    DataIterator dit = dbl.dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
      {
        phi[dit()].setVal(1.0);
        lph[dit()].setVal(0.0);
      }

    Real time = 0;
    int nstep = 0;
    while ((time < tfinal) && (nstep < nstepmax))
      {
        //advance the time and step counter
        time += dt;
        nstep++;
        cout << "nstep = " << nstep << "  time = " << time << endl;;
        //exchange ghost cell information
        phi.exchange(phi.interval());
        //do boundary conditions
        for (dit.reset(); dit.ok(); ++dit)
          {
            FArrayBox& soln = phi[dit()];
            const Box& region = dbl.get(dit());
            for (int idir = 0; idir < SpaceDim; idir++)
              {
                Box bcbox;
                if (region.smallEnd(idir) == domain.smallEnd(idir))
                  {
                    int ichop = domain.smallEnd(idir);
                    bcbox = soln.box();
                    bcbox.chop(idir, ichop);
                    int iside = -1;
                    FORT_BNDRYSUB(CHF_FRA1(soln,0),
                                  CHF_BOX(bcbox),
                                  CHF_INT(idir),
                                  CHF_INT(iside));
                  }
                if (region.bigEnd(idir) == domain.bigEnd(idir))
                  {
                    int ichop = domain.bigEnd(idir)+1;
                    Box chop_box = soln.box();
                    bcbox = chop_box.chop(idir,ichop);
                    int iside = 1;
                    FORT_BNDRYSUB(CHF_FRA1(soln,0),
                                  CHF_BOX(bcbox),
                                  CHF_INT(idir),
                                  CHF_INT(iside));
                  }
              }
          }
        //advance the solution in time
        for (dit.reset(); dit.ok(); ++dit)
          {
            FArrayBox& soln = phi[dit()];
            FArrayBox& lphi = lph[dit()];
            const Box& region = dbl.get(dit());
            FORT_HEATSUB(CHF_FRA1(soln,0),
                         CHF_FRA1(lphi,0),
                         CHF_BOX(region),
                         CHF_REAL(dt),
                         CHF_REAL(dx),
                         CHF_REAL(coeff));
          }
      } //end loop through time

    //output solution
#ifdef CH_USE_HDF5
    string filename("phi.hdf5");
    WriteUGHDF5(filename, dbl, phi, domain);
#else
    cout << "HDF5 is not defined, so I don't know how to write output!"
         << endl;
#endif

  }//end scoping trick
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return(0);
}
