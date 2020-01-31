#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// #include <cstdio>

//this one has getNearPeriodic
#include "PiecewiseLinearFillPatch.H"
#include "NWOQuadCFInterp.H"
#include "IntVectSet.H"
#include "DebugOut.H"
#include "NamespaceHeader.H"

void 
NWOQuadCFInterp::
define(/// layout at this level
       const DisjointBoxLayout&  a_thisDisjointBoxLayout,
       /// layout at coarser level
       const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
       /// number of variables
       const int&                a_numStates,
       /// problem domain on the coarser level
       const ProblemDomain&      a_coarseDomain,
       /// refinement ratio between this level and the coarser level
       const int&                a_refineCoarse,
       /// number of layers of ghost cells to fill by interpolation
       const int&                a_interpRadius)
{
  // Cache data
  m_numStates = a_numStates;
  m_coarseDomain = a_coarseDomain;
  m_refineCoarse = a_refineCoarse;
  m_interpRadius = a_interpRadius;

  m_layout = a_thisDisjointBoxLayout;
  m_coarseLayout = a_coarserDisjointBoxLayout;

  m_refineVect = m_refineCoarse * IntVect::Unit;

  ProblemDomain fineDomain = refine(m_coarseDomain, m_refineVect);

  coarsen(m_layoutCoarsened, m_layout, m_refineVect);

  int coarseGhost = a_interpRadius/a_refineCoarse  + 2;
  m_coarsenedFineData.define(m_layoutCoarsened, m_numStates, coarseGhost*IntVect::Unit);

  m_cfivs.define(m_layout);
  for(DataIterator dit = m_layout.dataIterator(); dit.ok(); ++dit)
    {
      Box ghostBox = m_layout[dit()];
      ghostBox.grow(a_interpRadius);
      ghostBox &= fineDomain;
      m_cfivs[dit].define(dit(), m_layout, ghostBox);
    }
  // Everything is defined now.
  m_defined = true;
}


//////////////////////////////////////////////////////////////////////////////
void 
NWOQuadCFInterp::
coarseFineInterp(/// interpolated solution on this level
                 LevelData<FArrayBox>&         a_fineData,
                 /// solution on coarser level
                 const LevelData<FArrayBox>&   a_coarseData,
                 /// starting coarse data component
                 int                           a_srcComp,
                 /// starting fine data component
                 int                           a_dstComp,
                 /// number of data components to interpolate
                 int                           a_numComp)
{
  CH_assert(m_defined);

  const Interval srcInterval(a_srcComp, a_srcComp + a_numComp-1);
  const Interval dstInterval(a_dstComp, a_dstComp + a_numComp-1);

  a_coarseData.copyTo(srcInterval,
                      m_coarsenedFineData,
                      dstInterval);

  int ibox = 0;
  for (DataIterator dit = m_layout.dataIterator();dit.ok(); ++dit)
    {
      interpOnPatch(a_fineData[dit()], m_coarsenedFineData[dit()], dit(), a_srcComp, a_dstComp, a_numComp);
      ibox++;
    }
}


//////////////////////////////////////////////////////////////////////////////
void 
NWOQuadCFInterp::
homogeneousCoarseFineInterp(/// interpolated solution on this level
                            LevelData<FArrayBox>&         a_fineData,
                            /// solution on coarser level
                            int                           a_srcComp,
                            /// starting fine data component
                            int                           a_dstComp,
                            /// number of data components to interpolate
                            int                           a_numComp)
{
  CH_assert(m_defined);

  for(DataIterator dit = m_coarsenedFineData.dataIterator(); dit.ok(); ++dit)
    {
      m_coarsenedFineData[dit()].setVal(0.);
    }

  for (DataIterator dit = m_layout.dataIterator();dit.ok(); ++dit)
    {
      interpOnPatch(a_fineData[dit()], m_coarsenedFineData[dit()], dit(), a_srcComp, a_dstComp, a_numComp);
    }
}


/////////////
void 
NWOQuadCFInterp::
interpOnPatch(FArrayBox&         a_fineData,
              const FArrayBox&   a_coarseData,
              const DataIndex&   a_dit,
              int                a_srcComp,
              int                a_dstComp,
              int                a_numComp)
{
  const IntVectSet& cfivsFine = m_cfivs[a_dit].getIVS();
  //dumpIVS(&cfivsFine);
  Real dxFine = 1;
  Real dxCoar = m_refineCoarse;
  for(IVSIterator ivsit(cfivsFine); ivsit.ok(); ++ivsit)
    {
      const IntVect& ivFine = ivsit();
      int ideb = 0;
      if((ivFine[0]==16) && (ivFine[1]==8))
        {
          ideb = 1;
        }
      IntVect ivCoar = coarsen(ivFine, m_refineCoarse);
      RealVect fineLoc, coarLoc;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          fineLoc[idir] = dxFine*(ivFine[idir] + 0.5);
          coarLoc[idir] = dxCoar*(ivCoar[idir] + 0.5);
        }

      for(int icomp = 0; icomp < a_numComp; icomp++)
        {
          int srcComp = a_srcComp + icomp;
          int dstComp = a_dstComp + icomp;

          RealVect firstDerivs, secondDerivs, mixedDerivs;
          getDerivs(firstDerivs, secondDerivs, mixedDerivs,
                    a_coarseData, ivCoar, dxCoar, icomp);

          Real coarValue = a_coarseData(ivCoar, srcComp);
          RealVect distance = fineLoc - coarLoc;
          Real fineValue;
          extrapolateValue(fineValue, coarValue, firstDerivs, secondDerivs, mixedDerivs,  distance);
          
          a_fineData(ivFine, dstComp) = fineValue;
        }
                           
    }
}
void 
NWOQuadCFInterp::
getDerivs(RealVect& firstDerivs, 
          RealVect& secondDerivs,
          RealVect& mixedDerivs, 
          const FArrayBox & a_data, 
          const IntVect& a_ivCoar,
          const Real   & a_dx,
          const int    & a_icomp)
{
  const IntVect& iv = a_ivCoar;
  int icomp = a_icomp;
  Real dx = a_dx;
  //single direction derivs--use centered diffs if possible, one-sided otherwise
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      IntVect ivhi = iv + BASISV(idir);
      IntVect ivlo = iv - BASISV(idir);
      bool hasHi = m_coarseDomain.contains(ivhi);
      bool hasLo = m_coarseDomain.contains(ivlo);
      if(hasHi && hasLo)
        {
          firstDerivs[idir]  = (a_data(ivhi, icomp) -  a_data(ivlo, icomp))/(2.*a_dx);
          secondDerivs[idir] = (a_data(ivhi, icomp) +  a_data(ivlo, icomp) - 2*a_data(iv, icomp))/(a_dx*a_dx);
        }
      else if(hasHi)
        {
          IntVect ivVeryHi = ivhi + BASISV(idir);
          firstDerivs[idir]  = (a_data(ivhi, icomp)     -  a_data(iv  , icomp))/(a_dx);
          secondDerivs[idir] = (a_data(ivVeryHi, icomp) +  a_data(iv  , icomp) - 2*a_data(ivhi, icomp))/(a_dx*a_dx);
        }
      else if(hasLo)
        {
          IntVect ivVeryLo = ivlo - BASISV(idir);
          firstDerivs[idir]  = (a_data(iv  , icomp)     -  a_data(ivlo, icomp))/(a_dx);
          secondDerivs[idir] = (a_data(ivVeryLo, icomp) +  a_data(iv  , icomp) - 2*a_data(ivlo, icomp))/(a_dx*a_dx);
        }
      else
        {
          firstDerivs[idir]  = 0;
          secondDerivs[idir] = 0;
        }
    }
  //now for that evil mixed deriv stuff --- 2d only has one, 3d has 3
  //this is to keep from doing each pair twice
  Vector<int> doneThisPair(SpaceDim, 0);
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      for(int jdir = 0; jdir < SpaceDim; jdir++)
        {
          if(idir != jdir)
            {
              int index = getMixedIndex(idir, jdir);
              if(doneThisPair[index] == 0)
                {
                  doneThisPair[index]  = 1;
                  IntVect ivhiI    = iv + BASISV(idir);
                  IntVect ivloI    = iv - BASISV(idir);
                  IntVect ivhiJ    = iv + BASISV(jdir);
                  IntVect ivloJ    = iv - BASISV(jdir);
                  IntVect ivhiIhiJ = iv + BASISV(idir) + BASISV(jdir);
                  IntVect ivloIloJ = iv - BASISV(idir) - BASISV(jdir);
                  IntVect ivloIhiJ = iv - BASISV(idir) + BASISV(jdir);
                  IntVect ivhiIloJ = iv + BASISV(idir) - BASISV(jdir);

                  bool hasIvhiIhiJ = m_coarseDomain.contains(ivhiIhiJ);
                  bool hasIvloIloJ = m_coarseDomain.contains(ivloIloJ);
                  bool hasIvloIhiJ = m_coarseDomain.contains(ivloIhiJ);
                  bool hasIvhiIloJ = m_coarseDomain.contains(ivhiIloJ);

                  //just go through the corners and compute each mixed deriv that you have
                  Real derivSum = 0;
                  int numDerivs = 0;
                  if(hasIvhiIhiJ)
                    {
                      Real dathiIhiJ = a_data(ivhiIhiJ, icomp);
                      Real dathiI    = a_data(ivhiI   , icomp);
                      Real dathiJ    = a_data(ivhiJ   , icomp);
                      Real datcen    = a_data(iv      , icomp);
                      Real mixedD = ((dathiIhiJ - dathiJ) - (dathiI   - datcen))/dx/dx;
                      derivSum += mixedD;
                      numDerivs++;
                    }
                  if(hasIvloIloJ)
                    {
                      Real datloIloJ = a_data(ivloIloJ, icomp);
                      Real datloI    = a_data(ivloI   , icomp);
                      Real datloJ    = a_data(ivloJ   , icomp);
                      Real datcen    = a_data(iv      , icomp);
                      Real mixedD = ((datcen - datloI) -  (datloJ   - datloIloJ))/dx/dx;
                      derivSum += mixedD;
                      numDerivs++;
                    }
                  if(hasIvhiIloJ)
                    {
                      Real dathiIloJ = a_data(ivhiIloJ, icomp);
                      Real dathiI    = a_data(ivhiI   , icomp);
                      Real datloJ    = a_data(ivloJ   , icomp);
                      Real datcen    = a_data(iv      , icomp);
                      Real mixedD = ((dathiI - datcen) - (dathiIloJ   - datloJ))/dx/dx;
                      derivSum += mixedD;
                      numDerivs++;
                    }

                  if(hasIvloIhiJ)
                    {
                      Real datloIhiJ = a_data(ivloIhiJ, icomp);
                      Real datloI    = a_data(ivloI   , icomp);
                      Real dathiJ    = a_data(ivhiJ   , icomp);
                      Real datcen    = a_data(iv      , icomp);
                      Real mixedD = ((dathiJ - datloIhiJ) - (datcen - datloI))/dx/dx;
                      derivSum += mixedD;
                      numDerivs++;
                    }
              
                  Real derivAvg = 0;
                  if(numDerivs > 0)
                    {
                      derivAvg = derivSum/numDerivs;
                    }
                  mixedDerivs[index]  = derivAvg;
                }
            }
        }
    }
}
void 
NWOQuadCFInterp::
extrapolateValue(Real           & a_fineValue, 
                 const Real     & a_coarValue, 
                 const RealVect & a_firstDerivs, 
                 const RealVect & a_secondDerivs, 
                 const RealVect & a_mixedDerivs, 
                 const RealVect & a_distance)
{
  a_fineValue = a_coarValue;
  //add in first and second derivative contributions
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      Real dx = a_distance[idir];
      a_fineValue +=          dx*a_firstDerivs[idir];
      a_fineValue += (dx*dx/2.)*a_secondDerivs[idir];
    }
  //now for the evil mixed derivatives
#if CH_SPACEDIM==2
  Real dxdy = a_distance[0]*a_distance[1];
  a_fineValue +=  dxdy*a_mixedDerivs[0];
#else
  for(int iindex = 0; iindex < SpaceDim; iindex++)
    {
      int idir, jdir;
      getMixedDerivDirections(idir, jdir, iindex);
      Real dx = a_distance[idir];
      Real dy = a_distance[jdir];
      a_fineValue +=  dx*dy*a_mixedDerivs[iindex];
    }
#endif
}



#include "NamespaceFooter.H"
