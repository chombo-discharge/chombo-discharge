#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MultiBlockCoordSys.H"
#include "BitSet.H"
#include "BoxIterator.H"
#include <cmath> // contains floor()

#include "NamespaceHeader.H"

MultiBlockCoordSys::MultiBlockCoordSys()
  :
  m_gotCoordSysVect(false),
  m_gotMappingBlocks(false),
  m_gotBoundaries(false),
  m_realDim(SpaceDim),
  m_gotBlockTransformations(false),
  m_useScriptN(true)
{
}

MultiBlockCoordSys::~MultiBlockCoordSys()
{
  if (m_gotBlockTransformations)
    {
      int nblocks = numBlocks();
      for (int blockSrc = 0; blockSrc < nblocks; blockSrc++)
        {
          delete m_validNum[blockSrc];
          delete m_validBlocks[blockSrc];
          delete m_validTransformations[blockSrc];
        }
    }
}

int MultiBlockCoordSys::whichBlock(const Box& a_box) const
{
  CH_assert(gotMappingBlocks()); // need m_mappingBlocks
  for (int n=0; n<numBlocks(); n++)
    {
      if (m_mappingBlocks[n].contains(a_box)) return n;
    }
  return -1;
}


int MultiBlockCoordSys::whichBlock(const IntVect& a_iv) const
{
  CH_assert(gotMappingBlocks()); // need m_mappingBlocks
  for (int n=0; n<numBlocks(); n++)
    {
      if (m_mappingBlocks[n].contains(a_iv)) return n;
    }
  return -1;
}


int MultiBlockCoordSys::whichBlockOverlap(const Box& a_box) const
{
  CH_assert(gotMappingBlocks()); // need m_mappingBlocks
  bool found = false;
  int blockNumber = -1;
  for (int n=0; n<numBlocks(); n++)
    {
      if (m_mappingBlocks[n].intersectsNotEmpty(a_box))
        {
          if (!found)
            { // first block that overlaps
              found = true;
              blockNumber = n;
            }
          else
            { // already found a different block
              blockNumber = -1;
            }
        }
    }
  return blockNumber;
}


int MultiBlockCoordSys::whichBlockBuffered(const IntVect& a_iv,
                                           int a_bufferSize) const
{
  CH_assert(gotMappingBlocks()); // need m_mappingBlocks
  for (int n=0; n<numBlocks(); n++)
    {
      Box expandedBlock = grow(m_mappingBlocks[n], a_bufferSize);
      if (expandedBlock.contains(a_iv)) return n;
    }
  return -1;
}


void MultiBlockCoordSys::initializeBlockTransformations()
{
  CH_assert(gotMappingBlocks()); // need m_mappingBlocks
  CH_assert(gotBoundaries()); // need m_boundaries
  // Define and fill in m_allBlockTransformations.
  int nblocks = numBlocks();
  m_allBlockTransformations.resize(nblocks);
  for (int blockSrc = 0; blockSrc < nblocks; blockSrc++)
    {
      Vector<IndicesTransformation>& destTransformations =
        m_allBlockTransformations[blockSrc];
      destTransformations.resize(nblocks);
      destTransformations.assign(IndicesTransformation::Undefined);

      // For each block:  do we have a path for it yet?
      BitSet destVisited(nblocks, false);
      // For each block:  have we followed all paths from it yet?
      BitSet destFollowed(nblocks, false);

      // blockSrc to itself
      destTransformations[blockSrc] = IndicesTransformation::Identity;
      destVisited.setTrue(blockSrc);

      for (int iter = 0; iter < nblocks; iter++)
        {
          if (!destVisited.isFull()) // some destinations still unvisited
            for (int blockVisited = 0; blockVisited < nblocks; blockVisited++)
              if (destVisited[blockVisited] && !destFollowed[blockVisited])
                {
                  // blockVisited has been visited, but not yet followed.
                  // So now follow all faces from blockVisited.
                  const IndicesTransformation& srcToVisited =
                    destTransformations[blockVisited];
                  int faceID = 0;
                  for (SideIterator sit; sit.ok(); sit.next())
                    {
                      // Side::LoHiSide side = sit();
                      for (int idir = 0; idir < SpaceDim; idir++)
                        {
                          const BlockBoundary& bb =
                            m_boundaries[blockVisited][faceID];
                          if (bb.isConformal() || bb.isMapped())
                            {
                              int blockNew = bb.neighbor();
                              if (!destVisited[blockNew])
                                { // Found a path to a new block!
                                  IndicesTransformation visitedToNew =
                                    bb.getTransformation();
                                  destTransformations[blockNew] =
                                    srcToVisited.compose(visitedToNew);
                                  destVisited.setTrue(blockNew);
                                }
                            }
                          faceID++;
                        }
                    }
                  // We have now followed all faces from blockVisited.
                  destFollowed.setTrue(blockVisited);
                }
        }
    }

  // Define and fill in m_validNum, m_validBlocks, m_validTransformations.

  // Set maxValid to SpaceDim! (factorial).
  int maxValid = 1;
  for (int i = 1; i <= SpaceDim; i++) maxValid *= i;

  m_validNum.resize(nblocks);
  m_validBlocks.resize(nblocks);
  m_validTransformations.resize(nblocks);
  Box offsetBox(-IntVect::Unit, IntVect::Unit);
  for (int blockSrc = 0; blockSrc < nblocks; blockSrc++)
    {
      m_validNum[blockSrc] =
        new BaseFab<int>(offsetBox, 1);
      m_validBlocks[blockSrc] =
        new BaseFab<int>(offsetBox, maxValid);
      m_validTransformations[blockSrc] =
        new BaseFab<IndicesTransformation>(offsetBox, maxValid);
      BaseFab<int>& validNum = *m_validNum[blockSrc];
      BaseFab<int>& validBlocks = *m_validBlocks[blockSrc];
      BaseFab<IndicesTransformation>& validTransformations =
        *m_validTransformations[blockSrc];
      for (BoxIterator bitOff(offsetBox); bitOff.ok(); ++bitOff)
        {
          IntVect offset = bitOff();
          // Set validBlocks(offset, :)
          // and validTransformations(offset, :)
          // and validNum(offset, 0), which says how many components
          // are in the other two.

          if (offset == IntVect::Zero)
            { // stay within central block
              validBlocks(offset, 0) = blockSrc;
              validTransformations(offset, 0) = IndicesTransformation::Identity;
              validNum(offset, 0) = 1;
            }
          else
            {
              // Set changedDirs to list of directions where there is change.
              Vector<int> changedDirs;
              // Set changedSigns to list of whether to change + or - in direction in changedDirs.
              Vector<int> changedSigns;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  if (offset[idir] != 0)
                    {
                      changedDirs.push_back(idir);
                      changedSigns.push_back(offset[idir]);
                    }
                }
              // numChanged >= 1 because offset != IntVect::Zero
              int numChanged = changedDirs.size();
              // allPaths contains all possible permutations of the
              // dimensions that change.
              // The number of paths is numPaths = numChanged!;
              // each path is of length numChanged.
              int numPaths = 1;
              for (int i = 1; i <= numChanged; i++) numPaths *= i;

              Vector< Vector<int> > allPathsDir; //(numPaths);
              Vector<int> vec(numChanged);
              for (int i = 0; i < numChanged; i++) vec[i] = i;
              // To allPathsDir, append all permutations of vec
              // (indexed into changedDirs),
              // beginning with vec[0].
              permute(vec, 0, numChanged, allPathsDir, changedDirs);
              CH_assert(allPathsDir.size() == numPaths);

              Vector< Vector<int> > allPathsSign(numPaths);
              for (int ipath = 0; ipath < numPaths; ipath++)
                {
                  for (int istep = 0; istep < numChanged; istep++)
                    {
                      int dirOfStep = allPathsDir[ipath][istep];
                      allPathsSign[ipath].push_back(offset[dirOfStep]);
                    }
                }

              Vector<int> validBlocksVect;
              Vector<IndicesTransformation> validTransformationsVect;
              for (int ipath = 0; ipath < numPaths; ipath++)
                {
                  const Vector<int>& thisPathDir = allPathsDir[ipath];
                  const Vector<int>& thisPathSign = allPathsSign[ipath];
                  // the path direction and sign as with the current block.
                  Vector<int> currentPathDir(thisPathDir);
                  Vector<int> currentPathSign(thisPathSign);
                  // the block you're at along the path.
                  int blockCurrent = blockSrc;
                  IndicesTransformation transformationCurrent =
                    IndicesTransformation::Identity;
                  bool allOK = true;
                  // Recall that numChanged is the length of thisPath.
                  for (int istep = 0; istep < numChanged; istep++)
                    {
                      int thisStepDir = currentPathDir[istep];
                      int thisStepSign = currentPathSign[istep];
                      CH_assert(abs(thisStepSign) == 1);
                      int faceID = (thisStepSign == -1) ?
                        thisStepDir : (thisStepDir + SpaceDim);
                      const BlockBoundary& bb =
                        m_boundaries[blockCurrent][faceID];
                      if (bb.isConformal() || bb.isMapped())
                        {
                          int blockNew = bb.neighbor();
                          const IndicesTransformation& transformationToNew =
                            bb.getTransformation();

                          blockCurrent = blockNew;
                          transformationCurrent =
                            transformationCurrent.compose(transformationToNew);

                          // Now update currentPathDir and currentPathSign
                          // for steps beyond istep.
                          for (int istepMore = istep+1; istepMore < numChanged; istepMore++)
                            {
                              IntVect oldStepVector =
                                currentPathSign[istepMore] *
                                BASISV(currentPathDir[istepMore]);
                              IntVect newStepVector =
                                transformationToNew.transformVectorFwd(oldStepVector);
                              for (int idirFind = 0; idirFind < SpaceDim; idirFind++)
                                {
                                  if (newStepVector[idirFind] != 0)
                                    {
                                      currentPathDir[istepMore] = idirFind;
                                      currentPathSign[istepMore] = newStepVector[idirFind];
                                    }
                                }
                            }
                        }
                      else
                        { // get out of here; cannot continue path
                          allOK = false;
                          break;
                        }
                    }
                  // Now the path is complete.
                  if (allOK)
                    {
                      validBlocksVect.push_back(blockCurrent);
                      validTransformationsVect.push_back(transformationCurrent);
                    }
                }
              // Now done with all paths.
              int numValid = validBlocksVect.size();
              validNum(offset, 0) = numValid;
              for (int ivalid = 0; ivalid < numValid; ivalid++)
                {
                  validBlocks(offset, ivalid) = validBlocksVect[ivalid];
                  validTransformations(offset, ivalid) =
                    validTransformationsVect[ivalid];
                }
            }
        }
    }

  m_gotBlockTransformations = true;
}

void
MultiBlockCoordSys::permute(Vector<int>& a_vec,
                            int a_start,
                            int a_length,
                            Vector< Vector<int> >& a_allPathsDir,
                            const Vector<int>& a_changedDirs)
{
  // To a_allPathsDir, append all permutations of a_vec
  // (indexed into a_changedDirs),
  // beginning with a_vec[a_start].
  if (a_start == a_length - 1)
    { // OK, a_vec is filled.
      Vector<int> thisPathDir(a_length);
      for (int i = 0; i < a_length; i++)
        {
          thisPathDir[i] = a_changedDirs[a_vec[i]];
        }
      a_allPathsDir.push_back(thisPathDir);
    }
  else
    { // a_start < a_length - 1
      for (int i = a_start; i < a_length; i++)
        {
          // Swap a_vec[i] with a_vec[a_start],
          // then call permute on a_start + 1.
          int tmp = a_vec[i];
          a_vec[i] = a_vec[a_start];
          a_vec[a_start] = tmp;
          permute(a_vec, a_start + 1, a_length, a_allPathsDir, a_changedDirs);
          // Undo swap of a_vec[i] with a_vec[a_start].
          a_vec[a_start] = a_vec[i];
          a_vec[i] = tmp;
        }
    }
}


RealVect
MultiBlockCoordSys::blockRemappingGeneral(int a_nDst,
                                          const RealVect& a_xiSrc,
                                          int a_nSrc) const
{
  RealVect xiDst;
  CH_assert(m_gotCoordSysVect);
  if (a_nDst == a_nSrc)
    {
      xiDst = a_xiSrc;
    }
  else
    {
      RealVect X = m_coordSysVect[a_nSrc]->realCoord(a_xiSrc);
      xiDst = m_coordSysVect[a_nDst]->mappedCoord(X);
    }
  return xiDst;
}

void
MultiBlockCoordSys::validBlocksCells(Vector<int>& a_validBlocks,
                                     Vector<IntVect>& a_validCells,
                                     int a_origBlock,
                                     const IntVect& a_origCell)
{
  const Box& origBlockBox = m_mappingBlocks[a_origBlock];
  if (origBlockBox.contains(a_origCell))
    {
      a_validBlocks.push_back(a_origBlock);
      a_validCells.push_back(a_origCell);
    }
  else
    {
      IntVect offset = IntVect::Zero;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (a_origCell[idir] < origBlockBox.smallEnd(idir))
            offset[idir] = -1;
          else if (a_origCell[idir] > origBlockBox.bigEnd(idir))
            offset[idir] = +1;
        }
      BaseFab<int>& validNumFab = *m_validNum[a_origBlock];
      BaseFab<int>& validBlocksFab = *m_validBlocks[a_origBlock];
      BaseFab<IndicesTransformation>& validTransformationsFab =
        *m_validTransformations[a_origBlock];
      for (int ind = 0; ind < validNumFab(offset, 0); ind++)
        {
          int validBlock = validBlocksFab(offset, ind);
          a_validBlocks.push_back(validBlock);

          const IndicesTransformation& transformation =
            validTransformationsFab(offset, ind);
          IntVect validIndices = transformation.transformFwd(a_origCell);
          a_validCells.push_back(validIndices);

          CH_assert(m_mappingBlocks[validBlock].contains(validIndices));
        }
    }
}

IntVect
MultiBlockCoordSys::cellAvoidingBoundaries(const IntVect&   a_cell,
                                           int              a_blockNum,
                                           int              a_radius) const
{
  CH_TIME("MultiBlockCoordSys::cellAvoidingBoundaries");
  IntVect avoidingCell = a_cell;
  const Tuple<BlockBoundary, 2*SpaceDim>& boundariesBlock =
    m_boundaries[a_blockNum];
  const Box& blockBox = m_mappingBlocks[a_blockNum];
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (boundariesBlock[idir].isDomainBoundary())
        { // Physical boundary at low end of this block.
          int minAvoiding = blockBox.smallEnd(idir) + a_radius;
          if (a_cell[idir] < minAvoiding)
            {
              avoidingCell[idir] = minAvoiding;
            }
        }
      if (boundariesBlock[SpaceDim + idir].isDomainBoundary())
        { // Physical boundary at high end of this block.
          int maxAvoiding = blockBox.bigEnd(idir) - a_radius;
          if (a_cell[idir] > maxAvoiding)
            {
              avoidingCell[idir] = maxAvoiding;
            }
        }
    }
  return avoidingCell;
}

void
MultiBlockCoordSys::keepInDomain(Box& a_bx,
                                 int a_blockNum) const
{
  const Tuple<BlockBoundary, 2*SpaceDim>& boundariesBlock =
    m_boundaries[a_blockNum];
  const Box& blockBox = m_mappingBlocks[a_blockNum];
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (boundariesBlock[idir].isDomainBoundary())
        { // Domain boundary here.  Remove everything below it.
          int edgeDir = blockBox.smallEnd(idir);
          if (a_bx.bigEnd(idir) < edgeDir)
            { // the whole a_bx is below the boundary;
              // return empty box
              a_bx = Box();
              return;
            }
          else if (a_bx.smallEnd(idir) < edgeDir)
            {
              a_bx.setSmall(idir, edgeDir);
            }
        }
      if (boundariesBlock[SpaceDim + idir].isDomainBoundary())
        { // Domain boundary here.  Remove everything above it.
          int edgeDir = blockBox.bigEnd(idir);
          if (a_bx.smallEnd(idir) > edgeDir)
            { // the whole a_bx is above the boundary;
              // return empty box
              a_bx = Box();
              return;
            }
          else if (a_bx.bigEnd(idir) > edgeDir)
            {
              a_bx.setBig(idir, edgeDir);
            }
        }
    }
}

ProblemDomain
MultiBlockCoordSys::problemDomain(const Box&  a_baseBox,
                                  int         a_numGhost) const
{
  int blockNumber = whichBlock(a_baseBox);
  Box bxGhosted = grow(a_baseBox, a_numGhost);
  keepInDomain(bxGhosted, blockNumber);
  ProblemDomain pd(bxGhosted); // non-periodic
  return pd;
}

ProblemDomain
MultiBlockCoordSys::blockDomainOfBox(const Box&  a_baseBox) const
{
  int blockNumber = whichBlock(a_baseBox);
  const Box& blockBox = m_mappingBlocks[blockNumber];
  ProblemDomain pd(blockBox); // non-periodic
  return pd;
}

VectorTransformation
MultiBlockCoordSys::vectorBlockTransformation(int a_nDst,
                                              const RealVect& a_xiSrc,
                                              int a_nSrc) const
{
  VectorTransformation tfmSrcToDst;
  if (a_nDst == a_nSrc)
    {
      tfmSrcToDst = VectorTransformation::Identity;
    }
  else
    {
      VectorTransformation tfmSrcMappedToReal =
        m_coordSysVect[a_nSrc]->vectorTransformationMappedToReal(a_xiSrc);
      RealVect xiDst = blockRemappingGeneral(a_nDst, a_xiSrc, a_nSrc);
      VectorTransformation tfmRealToDstMapped =
        m_coordSysVect[a_nDst]->vectorTransformationRealToMapped(xiDst);
      tfmSrcToDst = tfmSrcMappedToReal.compose(tfmRealToDstMapped);
    }
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
  return tfmSrcToDst;
}

VectorTransformation
MultiBlockCoordSys::vectorBlockTransformationCenter(int a_nDst,
                                                    const IntVect& a_ivSrc,
                                                    int a_nSrc) const
{
  RealVect xiSrc = m_coordSysVect[a_nSrc]->centerMappedCoordinates(a_ivSrc);
  return vectorBlockTransformation(a_nDst, xiSrc, a_nSrc);
}


void
MultiBlockCoordSys::vectorBlockTransformationCenter(BaseFab<VectorTransformation>& a_vtFab,
                                                    int a_nDst,
                                                    const Box& a_srcBox,
                                                    int a_nSrc) const
{
  for (BoxIterator bit(a_srcBox); bit.ok(); ++bit)
    {
      IntVect ivSrc = bit();
      a_vtFab(ivSrc, 0) =
        vectorBlockTransformationCenter(a_nDst, ivSrc, a_nSrc);
    }
}

VectorTransformation
MultiBlockCoordSys::vectorBlockTransformationAvg(int a_nDst,
                                                 const IntVect& a_ivSrc,
                                                 int a_nSrc) const
{
  VectorTransformation tfm =
    vectorBlockTransformationCenter(a_nDst, a_ivSrc, a_nSrc);
  Real frac = 1./24.;
  tfm *= 1.0 - (2 * SpaceDim) * frac;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (int ioff = -1; ioff <= 1; ioff += 2)
        {
          IntVect ivAdj = a_ivSrc + ioff * BASISV(idir);
          tfm += vectorBlockTransformationCenter(a_nDst, ivAdj, a_nSrc) * frac;
        }
    }
  return tfm;
}

void
MultiBlockCoordSys::vectorBlockTransformationAvg(BaseFab<VectorTransformation>& a_vtFab,
                                                  int a_nDst,
                                                  const Box& a_srcBox,
                                                  int a_nSrc,
                                                  Interval a_fixedDims) const
{
  CH_assert(a_vtFab.box().contains(a_srcBox));
  int nInterpDims = SpaceDim - a_fixedDims.size();
  Real frac = 1./24.;
  Real fracComp = 1.0 - (2 * nInterpDims) * frac;
  for (BoxIterator bit(a_srcBox); bit.ok(); ++bit)
    {
      IntVect ivSrc = bit();
      VectorTransformation tfm =
        vectorBlockTransformationCenter(a_nDst, ivSrc, a_nSrc);
      tfm *= fracComp;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if ( !a_fixedDims.contains(idir) )
            {
              VectorTransformation tfmLo =
                vectorBlockTransformationCenter(a_nDst,
                                                ivSrc - BASISV(idir),
                                                a_nSrc);
              VectorTransformation tfmHi =
                vectorBlockTransformationCenter(a_nDst,
                                                ivSrc + BASISV(idir),
                                                a_nSrc);
              tfm += (tfmLo + tfmHi) * frac;
            }
        }
      a_vtFab(ivSrc, 0) = tfm;
    }
}

Tuple<VectorTransformation, SpaceDim>
MultiBlockCoordSys::vectorBlockTransformationGrad(int a_nDst,
                                                  const IntVect& a_ivSrc,
                                                  int a_nSrc) const
{
  Tuple<VectorTransformation, SpaceDim> tfmGrad;
  // const RealVect& dxSrc = m_coordSysVect[a_nSrc]->dx();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      VectorTransformation tfmLo =
        vectorBlockTransformationCenter(a_nDst, a_ivSrc - BASISV(idir), a_nSrc);
      VectorTransformation tfmHi =
        vectorBlockTransformationCenter(a_nDst, a_ivSrc + BASISV(idir), a_nSrc);
      tfmGrad[idir] = (tfmHi - tfmLo) / 2.; // (2 * dxSrc[idir]);
    }
  return tfmGrad;
}

void
MultiBlockCoordSys::vectorBlockTransformationGrad(BaseFab<VectorTransformation>& a_vtGradFab,
                                                  int a_nDst,
                                                  const Box& a_srcBox,
                                                  int a_nSrc,
                                                  Interval a_fixedDims) const
{
  CH_assert(a_vtGradFab.box().contains(a_srcBox));
  for (BoxIterator bit(a_srcBox); bit.ok(); ++bit)
    {
      IntVect ivSrc = bit();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if ( !a_fixedDims.contains(idir) )
            {
              VectorTransformation tfmLo =
                vectorBlockTransformationCenter(a_nDst,
                                                ivSrc - BASISV(idir),
                                                a_nSrc);
              VectorTransformation tfmHi =
                vectorBlockTransformationCenter(a_nDst,
                                                ivSrc + BASISV(idir),
                                                a_nSrc);
              a_vtGradFab(ivSrc, idir) = (tfmHi - tfmLo) / 2.;
            }
        }
    }
}


void
MultiBlockCoordSys::vectorBlockAvgTransformation(BaseFab<VectorAvgTransformation>& a_vtAvgFab,
                                                 int a_nDst,
                                                 const Box& a_srcBox,
                                                 int a_nSrc,
                                                 Interval a_fixedDims) const
{
  CH_assert(a_vtAvgFab.box().contains(a_srcBox));
  int nInterpDims = SpaceDim - a_fixedDims.size();
  Real frac = 1./24.;
  Real fracComp = 1.0 - (2 * nInterpDims) * frac;
  for (BoxIterator bit(a_srcBox); bit.ok(); ++bit)
    {
      IntVect ivSrc = bit();
      VectorTransformation tfm =
        vectorBlockTransformationCenter(a_nDst, ivSrc, a_nSrc);
      tfm *= fracComp;
      Vector<VectorTransformation> tfmGrad;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if ( !a_fixedDims.contains(idir) )
            {
              VectorTransformation tfmLo =
                vectorBlockTransformationCenter(a_nDst,
                                                ivSrc - BASISV(idir),
                                                a_nSrc);
              VectorTransformation tfmHi =
                vectorBlockTransformationCenter(a_nDst,
                                                ivSrc + BASISV(idir),
                                                a_nSrc);
              tfm += (tfmLo + tfmHi) * frac;
              tfmGrad.push_back((tfmHi - tfmLo) / 2.);
            }
        }
      a_vtAvgFab(ivSrc, 0) = VectorAvgTransformation(tfm, tfmGrad);
    }
}


Vector<RealVect>
MultiBlockCoordSys::displacements(const Vector<RealVect>&   a_dstCoords,
                                  const Vector<int>&        a_dstBlocks,
                                  const RealVect&           a_srcCoords,
                                  int                       a_srcBlock) const
{
  int len = a_dstCoords.size();
  CH_assert(a_dstBlocks.size() == len);

  Vector<RealVect> returnVec(len);
  RealVect srcX = m_coordSysVect[a_srcBlock]->realCoord(a_srcCoords);
  for (int i = 0; i < len; i++)
    {
      const RealVect& dstXi = a_dstCoords[i];
      int dstBlock = a_dstBlocks[i];
      RealVect dstX = m_coordSysVect[dstBlock]->realCoord(dstXi);
      returnVec[i] = dstX - srcX;
    }
  return returnVec;
}

Real
MultiBlockCoordSys::distance(const RealVect&   a_srcCoords,
                             int               a_srcBlock,
                             const RealVect&   a_dstCoords,
                             int               a_dstBlock) const
{
  Vector<RealVect> dstCoordsVec(1, a_dstCoords);
  Vector<int> dstBlockVec(1, a_dstBlock);
  // Note that we may call displacements() function in derived class.
  Vector<RealVect> dispVec = displacements(dstCoordsVec, dstBlockVec,
                                           a_srcCoords, a_srcBlock);
  RealVect disp = dispVec[0];
  RealVect disp2 = disp * disp;
  Real dist = sqrt(disp2.sum());
  return dist;
}

void
MultiBlockCoordSys::separateVolFlux(LevelData<FluxBox>& a_flux) const
{
}

void
MultiBlockCoordSys::setAllBoundaries(BlockBoundary::btype a_btype)
{
  CH_assert(gotMappingBlocks());
  int nblocks = numBlocks();
  m_boundaries.resize(nblocks);
  for (int iblock = 0; iblock < nblocks; iblock++)
    {
      Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[iblock];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          blockBoundaries[idir].define(a_btype);
          blockBoundaries[idir + SpaceDim].define(a_btype);
        }
    }
}

void
MultiBlockCoordSys::setBoundaryFromFaces(int a_srcBlock,
                                         int a_srcFaceID,
                                         int a_dstBlock,
                                         int a_dstFaceID,
                                         IntVect a_sign)
{
  CH_assert(gotMappingBlocks()); // need m_mappingBlocks
  IndicesTransformation it;

  int srcDim = a_srcFaceID;
  Side::LoHiSide srcSide = Side::Lo;
  if (a_srcFaceID >= SpaceDim)
    {
      srcDim -= SpaceDim;
      srcSide = Side::Hi;
    }

  int dstDim = a_dstFaceID;
  Side::LoHiSide dstSide = Side::Lo;
  if (a_dstFaceID >= SpaceDim)
    {
      dstDim -= SpaceDim;
      dstSide = Side::Hi;
    }

  it.defineFromFaces(m_mappingBlocks[a_srcBlock], srcDim, srcSide,
                     m_mappingBlocks[a_dstBlock], dstDim, dstSide, a_sign);
  m_boundaries[a_srcBlock][a_srcFaceID].define(it, a_dstBlock);
}

bool
MultiBlockCoordSys::isDefined() const
{
  return (m_gotCoordSysVect &&
          m_gotMappingBlocks &&
          m_gotBoundaries &&
          m_gotBlockTransformations);
}

#include "NamespaceFooter.H"
