#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MultiCartesianCS.H"
#include "CartesianBlockCS.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

MultiCartesianCS::MultiCartesianCS()
{
  MayDay::Error("MultiCartesianCS constructor needs count and wraparound");
}


MultiCartesianCS::MultiCartesianCS(IntVect a_count,
                                   IntVect a_wraparound)
{
  m_count = a_count;

  m_wraparound = a_wraparound;

  m_blockIndicesBox.define(IntVect::Zero, m_count - IntVect::Unit);

  m_numBlocks = m_blockIndicesBox.numPts();

  // These will be filled in in define(), as will m_blockIndices.
  m_coordSysVect.resize(m_numBlocks, NULL);
  m_mappingBlocks.resize(m_numBlocks);
}


MultiCartesianCS::~MultiCartesianCS()
{
  if (m_gotCoordSysVect)
    {
      for (int iblock = 0; iblock < m_numBlocks; iblock++)
        {
          delete m_coordSysVect[iblock];
        }
    }
}


void
MultiCartesianCS::define(const ProblemDomain& a_levelDomain,
                         const RealVect& a_dx)
{
  const Box& levelBox = a_levelDomain.domainBox();
  CH_assert(levelBox.smallEnd() == IntVect::Zero);

  IntVect baseLengths = levelBox.size() / (2*m_count - IntVect::Unit);
  CH_assert((2*m_count - IntVect::Unit) * baseLengths == levelBox.size());
  int baseLength = baseLengths[0];
  CH_assert(baseLengths == baseLength * IntVect::Unit);

  IntVect baseLo = IntVect::Zero;
  IntVect baseHi = (baseLength - 1) * IntVect::Unit;
  Box baseBox(baseLo, baseHi);

  int blockIndex = 0;
  m_blockIndices.define(m_blockIndicesBox, 1);
  m_origin.resize(m_numBlocks);
  BoxIterator bit(m_blockIndicesBox);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      m_origin[blockIndex] = (2*baseLength) * iv;
      m_blockIndices(iv, 0) = blockIndex;
      blockIndex++;
    }
  CH_assert(blockIndex == m_numBlocks);

  for (int iblock = 0; iblock < m_numBlocks; iblock++)
    {
      m_mappingBlocks[iblock] = baseBox + m_origin[iblock];
    }
  m_gotMappingBlocks = true;

  /*
    Define block boundaries.
   */
  defineBoundaries();

  initializeBlockTransformations();

  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      int iblock = m_blockIndices(iv, 0);
      const Box& bx = m_mappingBlocks[iblock];
      m_coordSysVect[iblock] = new CartesianBlockCS(iblock, iv, a_dx, bx);
    }

  m_gotCoordSysVect = true;
}

void
MultiCartesianCS::defineBoundaries()
{
  CH_assert(gotMappingBlocks());
  // length in dimension 0 of base block (index 0)
  int baseLength = m_mappingBlocks[0].size(0);
  m_boundaries.resize(m_numBlocks);

  IndicesTransformation it;
  BoxIterator bit(m_blockIndicesBox);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      int iblock = m_blockIndices(iv, 0);
      Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[iblock];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          bool wrap = (m_wraparound[idir] == 1);
          int fullLength = (2*m_count[idir] - 1) * baseLength;

          // Off low end: define blockBoundaries[idir].
          BlockBoundary& blockBoundaryLo = blockBoundaries[idir];
          if (iv[idir] == 0)
            {
              if (wrap)
                {
                  // This dimension has range 0:m_count[idir]-1,
                  // so wrap around from 0 to m_count[idir]-1.
                  IntVect ivOther = iv;
                  ivOther[idir] = m_count[idir]-1;
                  int blockOther = m_blockIndices(ivOther, 0);
                  // Wrap around from low end of iblock to high end of blockOther.
                  it.defineFromTranslation(fullLength * BASISV(idir));
                  blockBoundaryLo.define(it, blockOther);
                }
              else // no wraparound in this dimension
                {
                  blockBoundaryLo = BlockBoundary::BOUNDARY;
                }
            }
          else
            {
              int blockOther = m_blockIndices(iv - BASISV(idir), 0);
              // Jump across from low end of iblock to high end of blockOther.
              it.defineFromTranslation(-baseLength * BASISV(idir));
              blockBoundaryLo.define(it, blockOther);
            }

          // Off high end: define blockBoundaries[idir + SpaceDim].
          BlockBoundary& blockBoundaryHi = blockBoundaries[idir + SpaceDim];
          if (iv[idir] == m_count[idir]-1)
            {
              if (wrap)
                {
                  // This dimension has range 0:m_count[idir]-1,
                  // so wrap around from m_count[idir]-1 to 0.
                  IntVect ivOther = iv;
                  ivOther[idir] = 0;
                  int blockOther = m_blockIndices(ivOther, 0);
                  // Wrap around from high end of iblock to low end of blockOther.
                  it.defineFromTranslation(-fullLength * BASISV(idir));
                  blockBoundaryHi.define(it, blockOther);
                }
              else // no wraparound in this dimension
                {
                  blockBoundaryHi = BlockBoundary::BOUNDARY;
                }
            }
          else
            {
              int blockOther = m_blockIndices(iv + BASISV(idir), 0);
              // Jump across from high end of iblock to low end of blockOther.
              it.defineFromTranslation(baseLength * BASISV(idir));
              blockBoundaryHi.define(it, blockOther);
            }
        }
    }

  m_gotBoundaries = true;
}


void
MultiCartesianCS::blockRemapping(RealVect& a_xi_valid,
                                 int& a_n_valid,
                                 const RealVect& a_xiSrc,
                                 int a_nSrc) const
{
  RealVect X = m_coordSysVect[a_nSrc]->realCoord(a_xiSrc);

  // X_valid in 0:m_count-1, so wrap around where possible.
  RealVect X_valid = X;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (m_wraparound[idir] == 1)
        {
          Real length = m_count[idir] * 1.;
          if (X_valid[idir] < 0.)
            {
              X_valid[idir] += length;
            }
          else if (X_valid[idir] > length)
            {
              X_valid[idir] -= length;
            }
        }
    }

  IntVect iv_valid;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      iv_valid[idir] = floor(X_valid[idir]);
      if (iv_valid[idir] < 0)
        {
          iv_valid[idir] = 0;
        }
      else if (iv_valid[idir] > m_count[idir]-1)
        {
          iv_valid[idir] = m_count[idir]-1;
        }
    }
  a_n_valid = m_blockIndices(iv_valid, 0);

  // a_xi_valid = m_coordSysVect[a_n_valid]->mappedCoord(X);
}

Vector<RealVect>
MultiCartesianCS::displacements(const Vector<RealVect>&   a_dstCoords,
                                const Vector<int>&        a_dstBlocks,
                                const RealVect&           a_srcCoords,
                                int                       a_srcBlock) const
{
  Vector<RealVect> disps =
    MultiBlockCoordSys::displacements(a_dstCoords, a_dstBlocks,
                                      a_srcCoords, a_srcBlock);
  int len = disps.size();
  for (int i = 0; i < len; i++)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (m_wraparound[idir] == 1)
            { // allow for wraparound in dimension idir
              Real length = m_count[idir] * 1.;
              if (disps[i][idir] > 1.)
                {
                  disps[i][idir] -= length;
                }
              else if (disps[i][idir] < -1.)
                {
                  disps[i][idir] += length;
                }
            }
        }
    }
  return disps;
}


MultiBlockCoordSys*
MultiCartesianCSFactory::getCoordSys(const ProblemDomain& a_levelDomain,
                                     const RealVect& a_dx) const
{
  MultiCartesianCS* coordSysPtr = new MultiCartesianCS(m_count, m_wraparound);
  coordSysPtr->define(a_levelDomain, a_dx);
  return (static_cast <MultiBlockCoordSys*> (coordSysPtr));
}

#include "NamespaceFooter.H"
