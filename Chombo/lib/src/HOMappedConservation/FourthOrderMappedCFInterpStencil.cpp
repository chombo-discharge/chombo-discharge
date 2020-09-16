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

#include "FourthOrderMappedCFInterpStencil.H"
#include "MappedCFInterpF_F.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
FourthOrderMappedCFInterpStencil::FourthOrderMappedCFInterpStencil()
{
  m_defined = false;
}

//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
FourthOrderMappedCFInterpStencil::~FourthOrderMappedCFInterpStencil()
{
}


//////////////////////////////////////////////////////////////////////////////
// Define the object
void FourthOrderMappedCFInterpStencil::define(
                                              const IntVect&    a_bdryOffset,
                                              int               a_degree,
                                              Interval          a_fixedDims)
{
  m_bdryOffset = a_bdryOffset;
  m_degree = a_degree;
  m_fixedDims = a_fixedDims;
  m_nInterpDims = SpaceDim - m_fixedDims.size();

  IntVect interpUnit = IntVect::Unit;
  for (int dirf = m_fixedDims.begin(); dirf <= m_fixedDims.end(); dirf++)
    {
      interpUnit[dirf] = 0;
    }

  { // Find m_stencilCells.
    m_stencilCells.clear();
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
        m_stencilCells.push_back(iv);
      }

    // OUTER SET:
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        if ( !m_fixedDims.contains(idir) )
          { // Add cells beyond innerBox, if possible, in dimension idir

            int bdryOffsetDir = m_bdryOffset[idir];
            if (bdryOffsetDir <= 0)
              { // point available on right
                IntVect iv = IntVect::Zero;
                iv[idir] = center[idir] + 2;
                m_stencilCells.push_back(iv);
              }
            if (bdryOffsetDir >= 0)
              { // point available on left
                IntVect iv = IntVect::Zero;
                iv[idir] = center[idir] - 2;
                m_stencilCells.push_back(iv);
              }
          }
      }
  }

  m_stencilSize = m_stencilCells.size();

  Box degreeBox(IntVect::Zero, m_degree * interpUnit);
  // number of nonzero powers:  sum > 0 and sum <= m_degree
  m_numNonzeroPowers = 0;
  for (BoxIterator bit(degreeBox); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      int sumPowers = iv.sum();
      if ( (sumPowers > 0) && (sumPowers <= m_degree) )
        { // IntVect iv passes the test, so count it.
          m_numNonzeroPowers++;
        }
    }

  m_indices.clear();
  m_indicesFull.clear();
  for (int icell = 0; icell < m_stencilSize; icell++)
    {
      IntVect stencilCell = m_stencilCells[icell];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          int coordinate = stencilCell[idir];
          m_indicesFull.push_back(coordinate);
          if ( !m_fixedDims.contains(idir))
            {
              m_indices.push_back(coordinate);
            }
        }
    }

  m_avgPowers.resize(m_stencilSize * m_numNonzeroPowers);
  FORT_STENCILAVGPOWERS(CHF_VR(m_avgPowers),
                        CHF_CONST_VI(m_indicesFull),
                        CHF_CONST_INT(m_stencilSize),
                        CHF_CONST_INT(m_numNonzeroPowers),
                        CHF_CONST_INT(m_degree),
                        CHF_BOX(degreeBox));

  // Everything is defined now.
  m_defined = true;
}

#include "NamespaceFooter.H"
