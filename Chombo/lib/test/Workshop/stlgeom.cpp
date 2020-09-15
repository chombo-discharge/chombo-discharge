#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "STLIF.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"

#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "DebugDump.H"
#include "EBLevelGrid.H"
#include "EBCellFactory.H"
#include "EBAMRIO.H"
#include "UsingNamespace.H"
                                                       // 


int makeGeometry(Box& a_domain,
                 Real& a_dx)
{
  int eekflag =  0;
  int maxbox = 1024;
  //parse input file

  int nx = 64;

  RealVect origin = -RealVect::Unit;
  RealVect probhi =  RealVect::Unit;
  IntVect lo =        IntVect::Zero;
  IntVect hi = (nx-1)*IntVect::Unit;

  a_domain.setSmall(lo);
  a_domain.setBig(hi);

  Real prob_hi = (probhi[0]- origin[0]);
  a_dx = prob_hi/nx;

  string stlfile("sphere.stl");

  STLIF implicit(stlfile.c_str(), STLIF::ASCII);


  GeometryShop workshop(implicit, 0, a_dx*RealVect::Unit);
  //this generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain, origin, a_dx, workshop, maxbox);

  return eekflag;
}

/****/
void  fillEBMF(LevelData<EBCellFAB>& a_mf, const EBLevelGrid& a_eblg, const Real& a_dx)
{
  for(DataIterator dit = a_mf.dataIterator(); dit.ok(); ++dit)
  {
    Box valid = a_eblg.getDBL()[dit()];
    EBISBox ebisBox = a_eblg.getEBISL()[dit()];
    IntVectSet ivsbox(valid);
    Real pi = 4.*atan(1.0);
    a_mf[dit()].setVal(0.0);
    for(VoFIterator vofit(ivsbox, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const IntVect & iv = vofit().gridIndex();
      Real val = 1.0;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        Real x = (iv[idir]+0.5)*a_dx;
        val *= sin(pi*x);
      }
      a_mf[dit()](vofit(), 0) = val;
    }
  }
}
/****/

int stlgeom()
{
  
  Box domain;
  Real dx;
  //make the initial geometry
  pout() << "making EBIS" << endl;
  makeGeometry(domain, dx);

  Vector<Box> boxes(1, domain);;
  Vector<int> procs(1, 0);
  DisjointBoxLayout dbl(boxes, procs);

  pout() << "filling a leveldata with data" << endl;
  EBLevelGrid eblg(dbl, domain, 0, Chombo_EBIS::instance());
  EBCellFactory fact(eblg.getEBISL());
  LevelData<EBCellFAB> mf(eblg.getDBL(), 1, IntVect::Zero, fact);
  fillEBMF(mf, eblg, dx);

  writeEBLevelName(mf, string("sphere.stl.hdf5"));
  return 0;
}


int
main(int argc,char **argv)
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int retval = 0;
  { 
#if CH_SPACEDIM==3

    retval = stlgeom();

  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#endif
  if(retval == 0)
  {
    pout() << "stl sphere test passed" << endl;
  }
  else
  {
    pout() << "stl sphere test failed with code " << retval << endl;
  }
  }
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return 0;
}


