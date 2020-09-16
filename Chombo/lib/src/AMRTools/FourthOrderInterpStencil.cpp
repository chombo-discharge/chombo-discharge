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

#include "FourthOrderInterpStencil.H"
#include "BoxIterator.H"
#include "LAPACKMatrix.H"

#include "NamespaceHeader.H"
//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
FourthOrderInterpStencil::FourthOrderInterpStencil()
{
  m_defined = false;
}

//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
FourthOrderInterpStencil::~FourthOrderInterpStencil()
{
}


//////////////////////////////////////////////////////////////////////////////
// Define the object
void FourthOrderInterpStencil::define(
                                      const IntVect&         a_bdryOffset,
                                      const int&             a_refineCoarse,
                                      const int&             a_degree,
                                      Interval               a_fixedDims)
{ CH_TIME("FourthOrderInterpStencil::define");
  m_bdryOffset = a_bdryOffset;
  m_refineCoarse = a_refineCoarse;
  m_degree = a_degree;
  m_fixedDims = a_fixedDims;

  IntVect interpUnit = IntVect::Unit;
  m_refineVect = m_refineCoarse * IntVect::Unit;
  for (int dirf = m_fixedDims.begin(); dirf <= m_fixedDims.end(); dirf++)
    {
      interpUnit[dirf] = 0;
      m_refineVect[dirf] = 1;
    }

  m_baseFineBox.define(IntVect::Zero, (m_refineCoarse - 1) * interpUnit);

  // index of zero cell in m_coarseBaseIndices; it will be in innerBox.
  int nbr0 = -1;

  { CH_TIME("finding stencil cells");
    // Find m_coarseBaseIndices.
    m_coarseBaseIndices.clear();
    // INNER SET:
    // First find the center of the inner set.
    IntVect center = IntVect::Zero;
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        if ( !m_fixedDims.contains(idir) )
          {
            switch (m_bdryOffset[idir])
              {
              case -1:
                { // boundary is immediately to left
                  center[idir] = 1;
                  break;
                }
              case 1:
                { // boundary is immediately to right
                  center[idir] = -1;
                  break;
                }
              }
          }
      }
    Box innerBox = Box(center - interpUnit,
                       center + interpUnit);
    for (BoxIterator bit(innerBox); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        if (iv == IntVect::Zero)
          {
            nbr0 = m_coarseBaseIndices.size();
          }
        m_coarseBaseIndices.push_back(iv);
      }

    // OUTER SET:
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        if ( !m_fixedDims.contains(idir) )
          { // Add cells beyond innerBox, if possible, in dimension idir

            int bdryOffsetDir = m_bdryOffset[idir];
            if (bdryOffsetDir <= 0)
              { // cell available on right
                IntVect iv = IntVect::Zero;
                iv[idir] = center[idir] + 2;
                m_coarseBaseIndices.push_back(iv);
              }
            if (bdryOffsetDir >= 0)
              { // cell available on left
                IntVect iv = IntVect::Zero;
                iv[idir] = center[idir] - 2;
                m_coarseBaseIndices.push_back(iv);
              }
          }
      }
  }
  m_stencilSize = m_coarseBaseIndices.size();
  CH_assert(nbr0 >= 0);
  CH_assert(nbr0 < m_stencilSize);
  CH_assert(m_coarseBaseIndices[nbr0] == IntVect::Zero);

  m_coarseToFineFab.define(m_baseFineBox, m_stencilSize);

  int numFinePoints = m_baseFineBox.numPts();

  Box degreeBox(IntVect::Zero, m_degree * interpUnit);

  Vector<IntVect> nonzeroPowersVect;
  for (BoxIterator bit(degreeBox); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      int sumPowers = iv.sum();
      if ( (sumPowers > 0) && (sumPowers <= m_degree) )
        { // IntVect iv passes the test, so include it.
          nonzeroPowersVect.push_back(iv);
        }
    }
  int nonzeroPowersNum = nonzeroPowersVect.size();

  // Fill in m_coarseToFineFab, which lives on
  // m_baseFineBox == (IntVect::Zero, (m_refineCoarse - 1) * interpUnit).

  LAPACKMatrix coarsePowers(m_stencilSize, nonzeroPowersNum);
  { // Fill coarsePowers.
    Real* coarsePowersData = coarsePowers.dataPtr();
    int indCoarsePow = 0;
    for (int ipow = 0; ipow < nonzeroPowersNum; ipow++)
      {
        const IntVect& powers = nonzeroPowersVect[ipow];
        for (int inbr = 0; inbr < m_stencilSize; inbr++)
          {
            const IntVect& nbr = m_coarseBaseIndices[inbr];
            coarsePowersData[indCoarsePow] = 1.;
            for (int idir = 0; idir < SpaceDim; idir++)
              {
                coarsePowersData[indCoarsePow] *=
                  power1dcoarseind0avg(nbr[idir], powers[idir]);
              }
            indCoarsePow++;
          }
      }
    CH_assert(indCoarsePow == m_stencilSize * nonzeroPowersNum);
  }

  LAPACKMatrix coarsePowersT(coarsePowers);
  coarsePowersT.transpose(); // dimensions (nonzeroPowersNum, m_stencilSize)

  LAPACKMatrix ATAinv;
  multiply(ATAinv, coarsePowersT, coarsePowers); // dimensions (nonzeroPowersNum, nonzeroPowersNum)
  ATAinv.invert();

  LAPACKMatrix coeffs;
  multiply(coeffs, ATAinv, coarsePowersT); // dimensions (nonzeroPowersNum, m_stencilSize);

  LAPACKMatrix finePowers(numFinePoints, nonzeroPowersNum);
  { // Fill finePowers.
    Real* finePowersData = finePowers.dataPtr();
    int indFinePow = 0;
    for (int ipow = 0; ipow < nonzeroPowersNum; ipow++)
      {
        const IntVect& powers = nonzeroPowersVect[ipow];
        for (BoxIterator bit(m_baseFineBox); bit.ok(); ++bit)
          {
            IntVect ivFine = bit();
            finePowersData[indFinePow] = 1.;
            for (int idir = 0; idir < SpaceDim; idir++)
              {
                finePowersData[indFinePow] *=
                  power1dfineind0avg(ivFine[idir], powers[idir]);
              }
            indFinePow++;
          }
      }
    CH_assert(indFinePow == numFinePoints * nonzeroPowersNum);
  }

  m_coarseToFineFab.define(m_baseFineBox, m_stencilSize);
  LAPACKMatrix cfInterpolator; // dimensions (numFinePoints, m_stencilSize)
  multiply(cfInterpolator, finePowers, coeffs);
  int ipoint = 0;
  for (BoxIterator bit(m_baseFineBox); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      for(int ivar = 0; ivar <  m_stencilSize; ivar++)
        {
          m_coarseToFineFab(iv, ivar) = cfInterpolator(ipoint, ivar);
        }
      ipoint++;
    }
  //  pout() << "cfInterpolator: ";
  //  cfInterpolator.pout();
  for (BoxIterator bit(m_baseFineBox); bit.ok(); ++bit)
    {
      IntVect ivFine = bit();
      Real sumCoeffs = 0.;
      for (int inbr = 0; inbr < m_stencilSize; inbr++)
        {
          sumCoeffs += m_coarseToFineFab(ivFine, inbr);
        }
      m_coarseToFineFab(ivFine, nbr0) += 1. - sumCoeffs;
    }

  // Everything is defined now.
  m_defined = true;
}

//////////////////////////////////////////////////////////////////////////////
void FourthOrderInterpStencil::fillFine(
                                        FArrayBox&         a_fineFab,
                                        const FArrayBox&   a_coarseFab,
                                        const IntVect&     a_coarseDataCell,
                                        const IntVect&     a_coarseToFineOffset) const
{
  CH_assert(m_defined);
  // This function fills the fine cells inside a_coarseDataCell
  // with interpolated data from a_coarseFab.
  // petermc, 3 Oct 2012, changed from m_refineCoarse to m_refineVect
  // because there is no refinement in m_fixedDims.
  IntVect fineBase = m_refineVect * a_coarseDataCell + a_coarseToFineOffset;
  // The fine cells to be filled in are at m_baseFineBox + fineBase,
  // where
  // m_baseFineBox == (IntVect::Zero, (m_refineCoarse - 1) * interpUnit).

  // The coarse cells from which we draw data are at
  // {a_coarseDataCell + m_coarseBaseIndices[i*D+[0:D-1]]: 0 <= i < m_stencilSize}.

  // added by petermc, 20 Aug 2009:  fill in a_fineFab only where you can.
  // (previously, had set shiftedFineBox = m_baseFineBox)
  Box shiftedFineBox(a_fineFab.box());
  shiftedFineBox.shift(-fineBase);
  shiftedFineBox &= m_baseFineBox; // intersect with [0:nref-1]^D or slice

  // Fill in a_fineFab(shiftedFineBox + fineBase)
  // using a_coarseFab(a_coarseDataCell + m_coarseBaseIndices)
  // and stencil coefficients in m_coarseToFineFab(shiftedFineBox).

  // Note that shiftedFineBox is contained in
  // m_baseFineBox == (IntVect::Zero, (m_refineCoarse - 1) * interpUnit).

  for (int icomp = 0; icomp < a_fineFab.nComp(); icomp++)
    {
      for (BoxIterator bit(shiftedFineBox); bit.ok(); ++bit)
        {
          IntVect ivFine = bit();
          Real val = 0.;
          for (int inbr = 0; inbr < m_stencilSize; inbr++)
            {
              IntVect ivCoarse = a_coarseDataCell + m_coarseBaseIndices[inbr];
              val +=
                m_coarseToFineFab(ivFine, inbr) * a_coarseFab(ivCoarse, icomp);
            }
          a_fineFab(fineBase + ivFine, icomp) = val;
        }
    }

  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused++;
}


//////////////////////////////////////////////////////////////////////////////
Real FourthOrderInterpStencil::power1d0avg(Real a_lower,
                                           Real a_upper,
                                           int a_ipwr) const
{
  Real avg = (pow(a_upper, a_ipwr+1) - pow(a_lower, a_ipwr+1)) /
    (a_upper - a_lower);
  if ((a_ipwr > 0) && (a_ipwr % 2 == 0))
    {
      avg -= pow(0.5, a_ipwr);
    }
  avg /= ((a_ipwr + 1) * 1.);
  return avg;
}


//////////////////////////////////////////////////////////////////////////////
Real FourthOrderInterpStencil::power1dcoarseind0avg(int a_offset,
                                                    int a_ipwr) const
{
  Real lower = a_offset*1. - 0.5;
  Real upper = a_offset*1. + 0.5;
  Real avg = power1d0avg(lower, upper, a_ipwr);
  return avg;
}


//////////////////////////////////////////////////////////////////////////////
Real FourthOrderInterpStencil::power1dfineind0avg(int a_offset,
                                                  int a_ipwr) const
{
  Real lower = a_offset*1. / Real(m_refineCoarse) - 0.5;
  Real upper = (a_offset + 1)*1. / Real(m_refineCoarse) - 0.5;
  Real avg = power1d0avg(lower, upper, a_ipwr);
  return avg;
}



#include "NamespaceFooter.H"
