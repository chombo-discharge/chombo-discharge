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

void outputSingleGrid(const FArrayBox& soln, const Box& domain);


int main(int argc, char* argv[])
{
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

  Box domain(IntVect::Zero, (nx-1)*IntVect::Unit);
  FArrayBox soln(grow(domain,1), 1);
  FArrayBox lphi(domain,1);

  //set phi to 1 for an initial condition
  soln.setVal(1.0);
  lphi.setVal(0.0);

  Real time = 0;
  int nstep = 0;
  while ((time < tfinal) && (nstep < nstepmax))
    {
      //advance the time and step counter
      time += dt;
      nstep++;
      cout << "nstep = " << nstep << "  time = " << time << endl;;
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
                   domain.loVect(), domain.hiVect(),
                   domain.loVect(), domain.hiVect(),
                   &dt, &dx, &coeff);
    } //end loop through time

  //output solution
  outputSingleGrid(soln, domain);
  return 0;
}

void outputSingleGrid(const FArrayBox& soln, const Box& domain)
{
  string filename("phi.hdf5");
  Vector<Box> vbox(1, domain);
  Vector<int> procID(1,0);
  DisjointBoxLayout dbl(vbox, procID);
  LevelData<FArrayBox> phi(dbl,1);
  DataIterator dit = dbl.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      phi[dit()].copy(soln);
    }
#ifdef CH_USE_HDF5
  WriteUGHDF5(filename, dbl, phi, domain);
#endif
}
