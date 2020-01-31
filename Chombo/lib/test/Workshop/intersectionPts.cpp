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

#include "SphereIF.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "DebugDump.H"
#include "EBLevelGrid.H"
#include "IFData.H"
#include "EBArith.H"
#include "UsingNamespace.H"


int makeGeometry(Box      & a_domain,
                 Real     & a_dx,
                 SphereIF*& a_spherePtr)
{
  int eekflag =  0;

  RealVect origin  = RealVect::Zero;
  int ncell = 8;
  Real domLen = 1;
  Real radius = 0.4;
  bool insideRegular = true;
  RealVect center = 0.5*RealVect::Unit;

  IntVect ivlo = IntVect::Zero;
  IntVect ivhi = (ncell-1)*IntVect::Unit;
  a_domain = Box(ivlo, ivhi);

  a_dx = domLen/ncell;

  a_spherePtr = new SphereIF(radius,center,insideRegular);

  RealVect vectDx = a_dx*RealVect::Unit;
  int verbosity = 0;
  GeometryShop workshop(*a_spherePtr,verbosity,vectDx);
  // this generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain, origin, a_dx, workshop);
  return eekflag;
}

int testIFDataIO()
{
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  RealVect origin;
  Real dx;
  Box domain;
  SphereIF* spherePtr = NULL;
  int eekflag = 0;
  eekflag = makeGeometry(domain,
                         dx,
                         spherePtr);

  int degreeP = 2;
  IndexTM<Real, SpaceDim> vectdx;
  vectdx.setAll(dx);

  Vector<Box> vbox(1, domain);
  Vector<int> proc(1, 0);
  DisjointBoxLayout dbl(vbox, proc);
  EBLevelGrid eblg(dbl, domain, 0, ebisPtr);
  for(DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
  {
    const Box& valid = dbl[dit()];
    const EBISBox& ebis = eblg.getEBISL()[dit()];
    IntVectSet ivs(valid);
    for(VoFIterator vofit(ivs, ebis.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      RealVect vofLocation = EBArith::getVoFLocation(vof, dx, RealVect::Zero);
      RvgDim cellCenter;
      for (int idir = 0;idir < SpaceDim; ++idir)
      {
        cellCenter[idir] = vofLocation[idir];
      }
      
      IFData<SpaceDim> ifdata(*spherePtr, vectdx, cellCenter, degreeP);
      ifdata.print(pout());
    }
  }

  
  delete spherePtr;
  ebisPtr->clear();
  return 0;
}
/***************/

int
main(int argc,char **argv)
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  int retval = testIFDataIO();
  if(retval != 0)
  {
    pout() << "ifdata io test failed with code" << retval << endl;
  }
  else
  {
    pout() << "ifdata io test passed" << endl;
  }

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
}




