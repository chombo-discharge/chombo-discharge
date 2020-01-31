#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SphereIF.H"
#include "MFIndexSpace.H"

#include "EBISLayout.H"
#include "ParmParse.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "DebugDump.H"
#include "UsingNamespace.H"

int readGeometryInfo(Box& a_domain,
                     Real& a_dx,
                     RealVect& a_origin,
                     RealVect& a_center,
                     Real& a_radius);

/***************/
// define a sphere EBIS.
/***************/
int makeGeometry(MFIndexSpace& mfIndexSpace,
                 const Box& domain,
                 const Real& dx,
                 const RealVect& origin,
                 const RealVect& center,
                 const Real& radius);

/***************/
/***************/
void dumpmemoryatexit();
/***************/
/***************/

char iter_str[80];

int
main(int argc,char **argv)
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // begin forever present scoping trick
  {

    const char* in_file = "sphere.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    RealVect center;
    Real radius;
    RealVect origin;
    Real dx;
    Box domain;
    int eekflag = 0;
    MFIndexSpace a, b;
    MFIndexSpace *mfIndexSpace, *oldMFIndexSpace;
    mfIndexSpace = &a;
    oldMFIndexSpace = &b;

    readGeometryInfo(domain,
                     dx,
                     origin,
                     center,
                     radius);

    int steps= 5;
    int step = 0;

    while (step < steps)
    {
      eekflag = makeGeometry(*mfIndexSpace,
                             domain,
                             dx,
                             origin,
                             center,
                             radius);
      if (eekflag != 0)
        {
          pout() << "non zero eek detected = " << eekflag << endl;
          // MayDay::Error("problem in makeGeometry");
        }

      step++;

      center[0]+= dx/3.0;
      center[1]+= dx/3.0;
      radius -= dx/10.0;
      {
        MFIndexSpace* swap = oldMFIndexSpace;
        oldMFIndexSpace = mfIndexSpace;
        mfIndexSpace = swap;
      }
    }

    pout() << "\n sphere test passed" << endl;
  } // end scoping trick

#ifdef CH_MPI
    MPI_Finalize();
#endif

    return 0;
}

int readGeometryInfo(Box& a_domain,
                     Real& a_dx,
                     RealVect& a_origin,
                     RealVect& a_center,
                     Real& a_radius)
{

  // parse input file
  ParmParse pp;

  Vector<int> n_cell(SpaceDim);
  pp.getarr("n_cell",n_cell,0,SpaceDim);

  CH_assert(n_cell.size() == SpaceDim);

  IntVect lo = IntVect::Zero;
  IntVect hi;

  for (int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if (n_cell[ivec] <= 0)
        {
          pout() << " bogus number of cells input = " << n_cell[ivec];
          return(-1);
        }
      hi[ivec] = n_cell[ivec] - 1;
    }

  a_domain.setSmall(lo);
  a_domain.setBig(hi);

  Vector<Real> prob_lo(SpaceDim, 1.0);
  Real prob_hi;

  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.get("prob_hi",prob_hi);

  a_dx = (prob_hi-prob_lo[0])/n_cell[0];
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_origin[idir] = prob_lo[idir];
    }

  pp.get("radius",a_radius);

  // ParmParse doesn't get RealVects, so work-around with Vector<Real>
  Vector<Real> vectorCenter;
  pp.getarr("center",vectorCenter,0,SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_center[idir] = vectorCenter[idir];
    }

  return 0;
} // end read from file

int makeGeometry(MFIndexSpace& mfIndexSpace,
                 const Box& a_domain,
                 const Real& a_dx,
                 const RealVect& a_origin,
                 const RealVect& a_center,
                 const Real& a_radius)
{

  int eekflag = 0;
  Vector<GeometryService*> geometry(2, NULL);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  SphereIF insideIF(a_radius, a_center, true);
  GeometryShop workshop0(insideIF,0,vectDx);
  workshop0.m_phase = 0;
  geometry[0] = &workshop0;

  SphereIF outsideIF(a_radius, a_center, false);
  GeometryShop workshop1(outsideIF,0,vectDx);
  workshop1.m_phase=1;
  geometry[1] = &workshop1;

  mfIndexSpace.define(a_domain,a_origin,a_dx,geometry);

  return eekflag;
}
