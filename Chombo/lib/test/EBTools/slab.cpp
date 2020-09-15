#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "UsingNamespace.H"

/***************/
/***************/
int makeEBISL(EBISLayout& a_ebisl,
              const DisjointBoxLayout& a_grids,
              const Box& a_domain,
              const int& a_nghost)
{
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
 CH_assert(ebisPtr->isDefined());
  ebisPtr->fillEBISLayout(a_ebisl, a_grids, a_domain, a_nghost);
  return 0;
}
/***************/
/***************/

int
makeLayout(DisjointBoxLayout& a_dbl,
           const Box& a_domain, int a_maxsize=0)
{
  ParmParse ppGeom("geom");

  int eekflag= 0;
  int ipieces;
  ipieces = Max(ipieces, 1);
  int maxsize = a_maxsize;
  if (maxsize == 0)
    ppGeom.get("maxboxsize",maxsize);

  Vector<Box> vbox(1, a_domain);
  domainSplit(a_domain, vbox,  maxsize);
  if (eekflag != 0)
    {
      pout() << "problem in domainsplit" << endl;
      return eekflag;
    }
  Vector<int>  procAssign;
  eekflag = LoadBalance(procAssign,vbox);
  if (eekflag != 0)
    {
      pout() << "problem in loadbalance" << endl;
      return eekflag;
    }
  a_dbl.define(vbox, procAssign);
  return eekflag;
}
/**********/
/**********/
int makeGeometry(Box& a_domain, Real& a_dx, Box& coveredBox)
{

  ParmParse ppGeom("geom");

  vector<int> n_cell(SpaceDim);
  ppGeom.getarr("n_cell",n_cell,0,SpaceDim);

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

  vector<Real> prob_lo(SpaceDim, 0.0);
  vector<Real> prob_hi(SpaceDim, 1.0);
  ppGeom.getarr("prob_lo",prob_lo,0,SpaceDim);
  ppGeom.getarr("prob_hi",prob_hi,0,SpaceDim);
  a_dx = (prob_hi[0]-prob_lo[0])/n_cell[0];

  vector<int> slab_lo(SpaceDim);
  ppGeom.getarr("slab_lo",slab_lo,0,SpaceDim);
  vector<int> slab_hi(SpaceDim);
  ppGeom.getarr("slab_hi",slab_hi,0,SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      lo[idir] = slab_lo[idir];
      hi[idir] = slab_hi[idir];
    }
  coveredBox = Box(lo,hi);
  SlabService slab(coveredBox);
  //this generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  RealVect origin = RealVect::Zero;
  ebisPtr->define(a_domain, origin, a_dx, slab, 8);

  return 0;
}
/***************/

