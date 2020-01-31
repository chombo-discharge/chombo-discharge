#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DoubleCartesianCS.H"
#include "CartesianBlockCS.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

DoubleCartesianCS::DoubleCartesianCS()
{
  m_coordSysVect.resize(NUMBLOCKS, NULL);

  m_mappingBlocks.resize(NUMBLOCKS);
}

DoubleCartesianCS::~DoubleCartesianCS()
{
  if (m_gotCoordSysVect)
    {
      for (int iblock = 0; iblock < NUMBLOCKS; iblock++)
        {
          delete m_coordSysVect[iblock];
        }
    }
}

void
DoubleCartesianCS::define(const ProblemDomain& a_levelDomain,
                          const RealVect& a_dx)
{
  const Box& levelBox = a_levelDomain.domainBox();
  int baseLength = levelBox.size(0) / 3;
  CH_assert(levelBox.size() == 3 * baseLength * IntVect::Unit);
  CH_assert(levelBox.smallEnd() == IntVect::Zero);

  IntVect baseLo = IntVect::Zero;
  IntVect baseHi = (baseLength - 1) * IntVect::Unit;
  Box baseBox(baseLo, baseHi);

  int blockIndex = 0;
  m_blockIndicesBox.define(IntVect::Zero, IntVect::Unit);
  m_blockIndices.define(m_blockIndicesBox, 1);
  BoxIterator bit(m_blockIndicesBox);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      m_origin[blockIndex] = (2*baseLength) * iv;
      m_blockIndices(iv, 0) = blockIndex;
      blockIndex++;
    }

  for (int iblock = 0; iblock < NUMBLOCKS; iblock++)
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
DoubleCartesianCS::defineBoundaries()
{
  CH_assert(gotMappingBlocks());
  // length in dimension 0 of central block (index 0)
  int baseLength = m_mappingBlocks[0].size(0);
  m_boundaries.resize(NUMBLOCKS);

  IndicesTransformation it;
  BoxIterator bit(m_blockIndicesBox);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      int iblock = m_blockIndices(iv, 0);
      Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[iblock];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          // Define blockBoundaries[idir] (off low end)
          // and blockBoundaries[idir + SpaceDim] (off high end).
          switch (iv[idir])
            {
            case 0:
              {
                int blockOther = m_blockIndices(iv + BASISV(idir), 0);
                // Wrap around from low end of iblock to high end of blockOther.
                it.defineFromTranslation(3 * baseLength * BASISV(idir));
                blockBoundaries[idir].define(it, blockOther);
                // Jump across from high end of iblock to low end of blockOther.
                it.defineFromTranslation(baseLength * BASISV(idir));
                blockBoundaries[idir + SpaceDim].define(it, blockOther);
              }
              break;
            case 1:
              {
                int blockOther = m_blockIndices(iv - BASISV(idir), 0);
                // Jump across from low end of iblock to high end of blockOther.
                it.defineFromTranslation(-baseLength * BASISV(idir));
                blockBoundaries[idir].define(it, blockOther);
                // Wrap around from high end of iblock to low end of blockOther.
                it.defineFromTranslation(-3*baseLength * BASISV(idir));
                blockBoundaries[idir + SpaceDim].define(it, blockOther);
              }
              break;
            default:
              MayDay::Error("DoubleCartesianCS::defineBoundaries bad index");
            }
        }
    }

  m_gotBoundaries = true;
}


void
DoubleCartesianCS::blockRemapping(RealVect& a_xi_valid,
                                  int& a_n_valid,
                                  const RealVect& a_xiSrc,
                                  int a_nSrc) const
{
  RealVect X = m_coordSysVect[a_nSrc]->realCoord(a_xiSrc);

  // X_valid in [0:2]^D.
  RealVect X_valid = X;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (X_valid[idir] < 0.)
        {
          X_valid[idir] += 2.;
        }
      else if (X_valid[idir] > 2.)
        {
          X_valid[idir] -= 2.;
        }
    }

  a_xi_valid = X_valid;
  IntVect iv_valid;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (X_valid[idir] > 1.0)
        {
          iv_valid[idir] = 1;
          a_xi_valid[idir] = X_valid[idir] + 1.;
        }
      else
        {
          iv_valid[idir] = 0;
          a_xi_valid[idir] = X_valid[idir];
        }
    }
  a_n_valid = m_blockIndices(iv_valid, 0);

  // a_xi_valid = m_coordSysVect[a_n_valid]->mappedCoord(X);
}

Vector<RealVect>
DoubleCartesianCS::displacements(const Vector<RealVect>&   a_dstCoords,
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
        { // allow for wraparound in dimension idir
          if (disps[i][idir] > 1.)
            disps[i][idir] -= 2.;
          else if (disps[i][idir] < -1.)
            disps[i][idir] += 2.;
        }
    }
  return disps;
}


void
DoubleCartesianCS::separateVolFlux(LevelData<FluxBox>& a_flux) const
{
  const DisjointBoxLayout& layout = a_flux.disjointBoxLayout();
  DataIterator dit = layout.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box bx = layout[dit];
      int blockNum = whichBlock(bx);
      IntVect blockOrigin = m_origin[blockNum];

      FluxBox& flub = a_flux[dit];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (blockOrigin[idir] > 0)
            { // not the first
              flub[idir].negate();
            }
        }
    }
}


MultiBlockCoordSys*
DoubleCartesianCSFactory::getCoordSys(const ProblemDomain& a_levelDomain,
                                      const RealVect& a_dx) const
{
  DoubleCartesianCS* coordSysPtr = new DoubleCartesianCS();
  coordSysPtr->define(a_levelDomain, a_dx);
  return (static_cast <MultiBlockCoordSys*> (coordSysPtr));
}

#include "NamespaceFooter.H"
