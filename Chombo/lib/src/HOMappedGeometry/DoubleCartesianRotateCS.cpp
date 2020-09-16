#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DoubleCartesianRotateCS.H"
#include "CartesianBlockRotateCS.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

DoubleCartesianRotateCS::DoubleCartesianRotateCS()
{
  m_coordSysVect.resize(NUMBLOCKS, NULL);

  m_mappingBlocks.resize(NUMBLOCKS);
}

DoubleCartesianRotateCS::~DoubleCartesianRotateCS()
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
DoubleCartesianRotateCS::define(const ProblemDomain& a_levelDomain,
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

  m_perm.define(m_blockIndicesBox, 1);
  m_sign.define(m_blockIndicesBox, 1);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      // int iblock = m_blockIndices(iv, 0);

      IntVect permBlock = IntVect(D_DECL6(0,1,2,3,4,5));
      IntVect signBlock = IntVect::Unit;
#if CH_SPACEDIM >= 2
      if (iv[0] != iv[1])
        { // swap dimensions
          permBlock[0] = 1; permBlock[1] = 0;
        }
      if (iv[0] == 1) signBlock[0] = -1;
      if (iv[1] == 1) signBlock[1] = -1;
      //      if (iv[0] == 1) // (iv == IntVect::Unit)
      //        {
      //          signBlock[0] = -1;
      //        }
#endif
      m_perm(iv, 0) = permBlock;
      m_sign(iv, 0) = signBlock;
    }

  m_vecTransform.define(m_blockIndicesBox, 1);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      int iblock = m_blockIndices(iv, 0);
      const Box& bx = m_mappingBlocks[iblock];

      IntVect permBlock = m_perm(iv, 0);
      IntVect signBlock = m_sign(iv, 0);
      CartesianBlockRotateCS*
        cBlockRotCS = new CartesianBlockRotateCS(iblock, iv, a_dx, bx,
                                                 permBlock, signBlock);
      m_coordSysVect[iblock] = cBlockRotCS;
      m_vecTransform(iv, 0) = cBlockRotCS->vecTransform();
    }

  /*
    Define block boundaries.
   */
  defineBoundaries();

  initializeBlockTransformations();

  m_gotCoordSysVect = true;
}

void
DoubleCartesianRotateCS::defineBoundaries()
{
  CH_assert(gotMappingBlocks());
  // length in dimension 0 of central block (index 0)
  // int baseLength = m_mappingBlocks[0].size(0);
  m_boundaries.resize(NUMBLOCKS);

  BoxIterator bit(m_blockIndicesBox);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect srcIV = bit();
      int srcBlock = m_blockIndices(srcIV, 0);
      Box srcBox = m_mappingBlocks[srcBlock];
      VectorTransformation srcVecTransform = m_vecTransform(srcIV, 0);
      int srcFaceID = 0;
      // CartesianBlockRotateCS* coordSysSrc =
      // dynamic_cast<CartesianBlockRotateCS*>(m_coordSysVect[srcBlock]);
      // CH_assert(coordSysSrc != NULL);
      for (SideIterator sit; sit.ok(); sit.next())
        {
          Side::LoHiSide side = sit();
          for (int srcDim = 0; srcDim < SpaceDim; srcDim++)
            {
              // unit vector pointing to face in mapped space
              RealVect srcPointing = RealVect(sign(side) * BASISV(srcDim));
              // unit vector pointing to face in rotated mapped space
              RealVect srcRotPointing = srcVecTransform.transform(srcPointing);
              // Get dimToOther and signToOther from srcRotPointing.

              int dimToOther = findDim(srcRotPointing);
              // int signToOther = findSign(srcRotPointing);
              // Find adjacent block.
              int upOrDown = (srcIV[dimToOther] == 0) ? 1 : -1;
              int dstBlock =
                m_blockIndices(srcIV + upOrDown*BASISV(dimToOther), 0);
              Box dstBox = m_mappingBlocks[dstBlock];

              // The unit vector srcRotPointing is the SAME unit vector
              // in the transformed destination box.
              // Need to rotate it back to dstPointing.
              CartesianBlockRotateCS* coordSysDst =
                dynamic_cast<CartesianBlockRotateCS*>(m_coordSysVect[dstBlock]);
              VectorTransformation dstVecTransform =
                coordSysDst->vecTransform();
              VectorTransformation dstInvVecTransform =
                dstVecTransform.inverse();
              RealVect dstPointing =
                dstInvVecTransform.transform(srcRotPointing);
              int dstDim = findDim(dstPointing);
              int dstSign = findSign(dstPointing);
              // dstPointing points FROM from dest face to center.
              // So dstSign == +1 means low face,
              // and dstSign == -1 means high face.
              int dstFaceID = dstDim;
              if (dstSign == -1) dstFaceID += SpaceDim;

              // Now figure out IntVect signBoundary.
              VectorTransformation compositeVecTransform =
                srcVecTransform.compose(dstInvVecTransform);
              IntVect signBoundary;
              for (int idir = 0; idir < SpaceDim; idir++)
                { // Set signBoundary[idir] to sign of compositeVecTransform
                  // on unit vector in direction idir.
                  RealVect evec = RealVect(BASISV(idir));
                  RealVect evecTransform =
                    compositeVecTransform.transform(evec);
                  // int evecDim = findDim(evecTransform);
                  int evecSign = findSign(evecTransform);
                  signBoundary[idir] = evecSign;
                }
              setBoundaryFromFaces(srcBlock, srcFaceID,
                                   dstBlock, dstFaceID, signBoundary);

              srcFaceID++;
            }
        }
    }

  m_gotBoundaries = true;
}


void
DoubleCartesianRotateCS::blockRemapping(RealVect& a_xi_valid,
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

  RealVect xi_valid_rot = X_valid;
  IntVect iv_valid;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (X_valid[idir] > 1.0)
        {
          iv_valid[idir] = 1;
          xi_valid_rot[idir] = X_valid[idir] + 1.;
        }
      else
        {
          iv_valid[idir] = 0;
          xi_valid_rot[idir] = X_valid[idir];
        }
    }
  a_n_valid = m_blockIndices(iv_valid, 0);

  // Now need to rotate back.
  CartesianBlockRotateCS* coordSysBlockValid =
    dynamic_cast<CartesianBlockRotateCS*>(m_coordSysVect[a_n_valid]);
  CH_assert(coordSysBlockValid != NULL);
  a_xi_valid = coordSysBlockValid->rotateBack(xi_valid_rot);
}

Vector<RealVect>
DoubleCartesianRotateCS::displacements(const Vector<RealVect>&   a_dstCoords,
                                       const Vector<int>&        a_dstBlocks,
                                       const RealVect&           a_srcCoords,
                                       int                       a_srcBlock) const
{
  // Calls m_coordSysVect[a_srcBlock]->realCoord(), which does rotations.
  Vector<RealVect> disps =
    MultiBlockCoordSys::displacements(a_dstCoords, a_dstBlocks,
                                      a_srcCoords, a_srcBlock);
  // In each dimension, use the wraparound displacement if it is shorter.
  // Should be in the range [-1:+1].
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
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
  return disps;
}


void
DoubleCartesianRotateCS::separateVolFlux(LevelData<FluxBox>& a_flux) const
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
        { // FIXME: Perhaps I need to use m_perm, m_sign?
          if (blockOrigin[idir] > 0)
            { // not the first
              flub[idir].negate();
            }
        }
    }
}

int
DoubleCartesianRotateCS::findDim(const RealVect& a_rv) const
{
  int retval = -1;
  bool found = false;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Real rvDir = a_rv[idir];
      if (abs(rvDir) > 0.5)
        {
          CH_assert(!found); // should not be found yet
          retval = idir;
          found = true;
        }
    }
  CH_assert(found);
  return retval;
}


int
DoubleCartesianRotateCS::findSign(const RealVect& a_rv) const
{
  int retval = 0;
  bool found = false;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Real rvDir = a_rv[idir];
      if (abs(rvDir) > 0.5)
        {
          CH_assert(!found); // should not be found yet
          retval = (rvDir > 0) ? 1 : -1;
          found = true;
        }
    }
  CH_assert(found);
  return retval;
}


MultiBlockCoordSys*
DoubleCartesianRotateCSFactory::getCoordSys(const ProblemDomain& a_levelDomain,
                                            const RealVect& a_dx) const
{
  DoubleCartesianRotateCS* coordSysPtr = new DoubleCartesianRotateCS();
  coordSysPtr->define(a_levelDomain, a_dx);
  return (static_cast <MultiBlockCoordSys*> (coordSysPtr));
}

#include "NamespaceFooter.H"
