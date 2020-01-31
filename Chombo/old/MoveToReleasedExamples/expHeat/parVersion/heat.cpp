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
#include "BoxIterator.H"
#include "Misc.H"
#include "Vector.H"
#include "REAL.H"
#include "Box.H"

#include "UsingNamespace.H"

#if CH_SPACEDIM==2

#ifdef CH_FORT_NOUNDERSCORE
#define FORT_HEATSUB heatsub2d
#else
#define FORT_HEATSUB heatsub2d_
#endif

extern "C"
{
  void FORT_HEATSUB(const Real* const phi,
                    const int*  const philo0,
                    const int*  const phihi0,
                    const int*  const philo1,
                    const int*  const phihi1,
                    const Real* const lph,
                    const int*  const lphlo0,
                    const int*  const lphhi0,
                    const int*  const lphlo1,
                    const int*  const lphhi1,
                    const int*  const boxlo,
                    const int*  const boxhi,
                    const int*  const domlo,
                    const int*  const domhi,
                    const Real* const dt,
                    const Real* const dx,
                    const Real* const nu);

}

#elif CH_SPACEDIM==3

#ifdef CH_FORT_NOUNDERSCORE
#define FORT_HEATSUB heatsub3d
#else
#define FORT_HEATSUB heatsub3d_
#endif


extern "C"
{
  void FORT_HEATSUB(const Real* const phi,
                    const int*  const philo0,
                    const int*  const phihi0,
                    const int*  const philo1,
                    const int*  const phihi1,
                    const int*  const philo2,
                    const int*  const phihi2,
                    const Real* const lph,
                    const int*  const lphlo0,
                    const int*  const lphhi0,
                    const int*  const lphlo1,
                    const int*  const lphhi1,
                    const int*  const lphlo2,
                    const int*  const lphhi2,
                    const int*  const boxlo,
                    const int*  const boxhi,
                    const int*  const domlo,
                    const int*  const domhi,
                    const Real* const dt,
                    const Real* const dx,
                    const Real* const nu);

}
#endif

void
makeGrids(Box& domain, DisjointBoxLayout& dbl, int& nx);

int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  {//scoping trick

    // number of points in each direction, domain length
    //and diffusion coefficient
    // make domain and grid layout
    Box domain;
    DisjointBoxLayout dbl;
    int nx;
    makeGrids(domain, dbl, nx);
    Real domainLen = 1.0;
    Real coeff = 1.0e-3;
    //define stopping conditions
    //set grid spacing and time step
    Real tfinal = 3.33;
    int nstepmax = 100;
    Real dx = domainLen/nx;
    Real dt = 0.8*dx*dx/(2.*SpaceDim*coeff);

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
        //advance the solution in time
        for (dit.reset(); dit.ok(); ++dit)
          {
            FArrayBox& soln = phi[dit()];
            FArrayBox& lphi = lph[dit()];
            const Box& region = dbl.get(dit());
            FORT_HEATSUB(soln.dataPtr(0),
                         &(soln.loVect()[0]), &(soln.hiVect()[0]),
                         &(soln.loVect()[1]), &(soln.hiVect()[1]),
#if CH_SPACEDIM==3
                         &(soln.loVect()[2]), &(soln.hiVect()[2]),
#endif
                         lphi.dataPtr(0),
                         &(lphi.loVect()[0]), &(lphi.hiVect()[0]),
                         &(lphi.loVect()[1]), &(lphi.hiVect()[1]),
#if CH_SPACEDIM==3
                         &(lphi.loVect()[2]), &(lphi.hiVect()[2]),
#endif
                         region.loVect(), region.hiVect(),
                         region.loVect(), region.hiVect(),
                         &dt, &dx, &coeff);
          }
      } //end loop through time

    //output solution
    string filename("phi.hdf5");
#ifdef CH_USE_HDF5
    WriteUGHDF5(filename, dbl, phi, domain);
#endif // CH_USE_HDF5

  }//end scoping trick
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return(0);
}
void
makeGrids(Box& domain, DisjointBoxLayout& dbl, int& len)
{
  int p = 4; // number of local patches in each direction.
  int Nx = 16; // number of grid points in each direction per local patch.

  len = Nx*p;
  Box    patch(IntVect::Zero, (Nx-1)*IntVect::Unit);
  Box skeleton(IntVect::Zero, (p-1)*IntVect::Unit);
  domain = Box(IntVect::Zero, (len-1)*IntVect::Unit);
  BoxIterator bit(skeleton);
  int thisProc = 0;
  Vector<Box> vbox;
  Vector<int> procId;
  for (bit.begin();bit.ok();++bit)
  {
    Box thisBox = patch + bit()*Nx;
    vbox.push_back(thisBox);
    procId.push_back(thisProc);
    thisProc = (thisProc + 1) % numProc();
  }

  dbl.define(vbox, procId);
  dbl.close();

  if (thisProc != 0)
  {
    cout << "warning: load is imbalanced, with " << thisProc <<
      "processors having one extra Box" << endl;
  }
}
