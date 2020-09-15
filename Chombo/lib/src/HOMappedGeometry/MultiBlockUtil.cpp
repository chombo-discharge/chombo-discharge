#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MultiBlockUtil.H"
#include "BoxIterator.H"
#include "MBStencilIterator.H"
#include "MBVectorStencilIterator.H"
#include "BoxCollapser.H"
#include "BoxFixedOff.H"
#include "MBMiscUtil.H"
#include "FourthOrderUtil.H"
#include "LAPACKMatrix.H"

#include "NamespaceHeader.H"

/// constructor
MultiBlockUtil::MultiBlockUtil()
{
  m_isDefined = false;
}

/// destructor
MultiBlockUtil::~MultiBlockUtil()
{
  undefine();
}

/// full constructor
MultiBlockUtil::MultiBlockUtil(const MultiBlockCoordSys*  a_coordSysPtr,
                               Interval                   a_fixedDims,
                               Vector<int>                a_fixedPt)
{
  define(a_coordSysPtr, a_fixedDims, a_fixedPt);
}

void 
MultiBlockUtil::undefine()
{
  if (m_isDefined)
    {
      m_interpDimsVect.clear();
      m_fixedDimsVect.clear();
      for (BoxIterator bitOff(m_offsetAllBox); bitOff.ok(); ++bitOff)
        {
          IntVect ivOff = bitOff();
          for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
            {
              delete m_transformsAll(ivOff, srcBlock);
            }
        }
    }
  m_isDefined = false;
}

void 
MultiBlockUtil::define(const MultiBlockCoordSys*  a_coordSysPtr,
                       Interval                   a_fixedDims,
                       Vector<int>                a_fixedPt)
{ CH_TIME("MultiBlockUtil::define");
  undefine();
  m_coordSysPtr = (MultiBlockCoordSys*) a_coordSysPtr;

  m_mappingBlocks = m_coordSysPtr->mappingBlocks();
  m_nblocks = m_mappingBlocks.size();

#if 0
  m_boundaries = m_coordSysPtr->boundaries();
#else
  // workaround until we figure out why the above breaks COGENT
  const Vector < Tuple< BlockBoundary, 2*SpaceDim > >& src_boundaries =
    m_coordSysPtr->boundaries();
  m_boundaries.resize(m_nblocks);
  for (int i=0; i<m_nblocks; ++i)
    {
      m_boundaries[i] = src_boundaries[i];
    }
#endif

  m_fixedDims = a_fixedDims;
  m_interpDimsVect.clear();
  m_fixedDimsVect.clear();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if ( a_fixedDims.contains(idir) )
        { // fix dimension idir
          m_interpUnit[idir] = 0;
          m_fixedDimsVect.push_back(idir);
        }
      else
        { // interpolate in dimension idir
          m_interpUnit[idir] = 1;
          m_interpDimsVect.push_back(idir);
        }
    }
  m_fixedPt = Vector<int>(a_fixedPt);

  // Use these in displacementPowers(), in cell-averaged case.
  int nInterpDims = m_interpDimsVect.size();
  m_convSize = 1+2*nInterpDims;
  m_convBaseCells.resize(m_convSize);
  m_convBaseCells.assign(IntVect::Zero);
  // Shifts are [0,0], [-1,0], [0,-1], [+1,0], [0,+1].
  // WAS for (int idir = 0; idir < SpaceDim; idir++)
  for (int ind = 0; ind < nInterpDims; ind++)
    {
      int idir = m_interpDimsVect[ind];
      m_convBaseCells[1 + ind][idir] -= 1;
      m_convBaseCells[1 + nInterpDims + ind][idir] += 1;
    }
  Real convScaling = 1./24.;
  m_convWeight.resize(m_convSize);
  m_convWeight.assign(convScaling);
  m_convWeight[0] = 1. - (2*nInterpDims)*convScaling;

  // Box m_offsetAllBox;
  m_offsetAllBox = Box(-IntVect::Unit, IntVect::Unit);

  //  m_transformsAll.resize(m_nblocks);
  m_transformsAll.define(m_offsetAllBox, m_nblocks);

  for (BoxIterator bitOff(m_offsetAllBox); bitOff.ok(); ++bitOff)
    {
      IntVect ivOff = bitOff();
      for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
        {
          m_transformsAll(ivOff, srcBlock) =
            new Vector<IndicesTransformation>(m_nblocks);
        }
    }

  for (int dstBlock = 0; dstBlock < m_nblocks; dstBlock++)
    {
      /*
        For each offset IntVect ivOff,
        validBlocks(ivOff, 0:validNum(ivOff, 0)-1)
        will hold the block numbers of valid blocks containing
        ghost cells in direction ivOff from block dstBlock,
        and validTransformations(ivOff, 0:validNum(ivOff, 0)-1)
        will hold the transformations from dstBlock to valid blocks
        in the same order.
      */
      const BaseFab<int>& validBlocks =
        m_coordSysPtr->validBlocks(dstBlock);
      const BaseFab<IndicesTransformation>& validTransformations =
        m_coordSysPtr->validTransformations(dstBlock);
      // validNum gives number of components of
      // validBlocks and validTransformations.
      const BaseFab<int>& validNum =
        m_coordSysPtr->validNum(dstBlock);
      for (BoxIterator bitOff(m_offsetAllBox); bitOff.ok(); ++bitOff)
        {
          IntVect ivOff = bitOff();
          int validNumOff = validNum(ivOff, 0);
          for (int i = 0; i < validNumOff; i++)
            {
              int srcBlock = validBlocks(ivOff, i);

              // Transformation from dstBlock to srcBlock.
              IndicesTransformation tfm = validTransformations(ivOff, i);
              // (*m_transformsAll[srcBlock])(ivOff, dstBlock) = tfm;
              (*m_transformsAll(ivOff, srcBlock))[dstBlock] = tfm;
            }
        }
    }

  m_isDefined = true;
}


IntVect
MultiBlockUtil::whichOffset(int             a_block,
                            const IntVect&  a_cell) const
{
  CH_assert(m_isDefined);
  const Box& blockBox = m_mappingBlocks[a_block];
  IntVect offset = IntVect::Zero;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (a_cell[idir] < blockBox.smallEnd(idir))
        offset[idir] = -1;
      else if (a_cell[idir] > blockBox.bigEnd(idir))
        offset[idir] = +1;
    }
  return offset;
}


const IndicesTransformation&
MultiBlockUtil::blockTransform(int a_srcBlock,
                               int a_dstBlock,
                               const IntVect& a_offset) const
{
  CH_assert(m_isDefined);
  const Vector<IndicesTransformation>& tfmSrcOff =
    transformsSrcOff(a_srcBlock, a_offset);
  return tfmSrcOff[a_dstBlock];
}


const Vector<IndicesTransformation>&
MultiBlockUtil::transformsSrcOff(int a_srcBlock,
                                 const IntVect& a_offset) const
{
  CH_assert(m_isDefined);
  const Vector<IndicesTransformation>& tfmSrcOff =
    *m_transformsAll(a_offset, a_srcBlock);
  return tfmSrcOff;
}


void
MultiBlockUtil::growIVS(IntVectSet& a_ivs,
                        int a_bufferSize) const
{
  CH_assert(m_isDefined);
  CH_TIME("MultiBlockUtil::growIVS");
  IntVectSet ivsGrown(a_ivs);
  ivsGrown.grow(a_bufferSize);
  for (IVSIterator ivsit(ivsGrown); ivsit.ok(); ++ivsit)
    {
      IntVect iv = ivsit();
      if ( !a_ivs.contains(iv) )
        {
          int blockValid = m_coordSysPtr->whichBlock(iv);
          if (blockValid != -1)
            { // iv is a valid cell of blockValid
              a_ivs |= iv;
            }
          else
            { // iv is not a valid cell of any block
              int dstBlock = m_coordSysPtr->whichBlockBuffered(iv, a_bufferSize);
              IntVect offset = whichOffset(dstBlock, iv);
              CH_assert(offset != IntVect::Zero);
              for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
                {
                  // tfm transforms from dstBlock to srcBlock,
                  // taking offset from dstBlock.
                  const IndicesTransformation& tfm = blockTransform(srcBlock,
                                                                    dstBlock,
                                                                    offset);
                  if ( tfm.isDefined() )
                    {
                      IntVect ivValid = tfm.transform(iv);
                      a_ivs |= ivValid;
                    }
                }
            }
        }
    }
}

void
MultiBlockUtil::commonVertexCells(Vector<IntVect>& a_neighborCells,
                                  Vector<int>& a_neighborBlocks,
                                  const IntVect& a_baseCell,
                                  int a_baseBlock) const
{
  CH_assert(m_isDefined);
  CH_TIME("MultiBlockUtil::commonVertexCells");
  Box baseBlockBox = m_mappingBlocks[a_baseBlock];

  // Check for the easy case: a_baseCell is away from all block boundaries.
  Box innerBlockBox = grow(baseBlockBox, -1);
  if (innerBlockBox.contains(a_baseCell))
    {
      a_neighborCells.clear();
      a_neighborBlocks.clear();
      Box neighborsBox(a_baseCell - IntVect::Unit,
                       a_baseCell + IntVect::Unit);
      for (BoxIterator bitNbr(neighborsBox); bitNbr.ok(); ++bitNbr)
        {
          IntVect nbr = bitNbr();
          a_neighborCells.push_back(nbr);
          a_neighborBlocks.push_back(a_baseBlock);
        }
    }
  else
    {
      Box offsetBox(IntVect::Zero, IntVect::Unit);
      Vector<IntVect> allVertices;
      Vector<int> allVertexBlocks;
      for (BoxIterator bitOff(offsetBox); bitOff.ok(); ++bitOff)
        {
          IntVect offset = bitOff();
          IntVect vertex = a_baseCell + offset;
          // Note that as a NODE, vertex is still valid within a_baseBlock.
          // First, find all aliases of vertex.
          Vector<IntVect> aliasVertices;
          Vector<int> aliasVertexBlocks;
          vertexAliases(aliasVertices, aliasVertexBlocks,
                        vertex, a_baseBlock);
          allVertices.append(aliasVertices);
          allVertexBlocks.append(aliasVertexBlocks);
        }
      validCellsFromVertices(a_neighborCells, a_neighborBlocks,
                             allVertices, allVertexBlocks);
    }
  return;
}

void
MultiBlockUtil::commonVertexStencilElements(Vector<MBStencilElement>& a_neighbors,
                                            const MBStencilElement& a_base) const
{
  CH_assert(m_isDefined);
  CH_TIME("MultiBlockUtil::commonVertexCells");
  int baseBlock = a_base.block();
  const IntVect& baseCell = a_base.cell();
  Box baseBlockBox = m_mappingBlocks[baseBlock];

  Box innerBlockBox = grow(baseBlockBox, -1);
  if (innerBlockBox.contains(baseCell))
    { // easy case: baseCell is away from all block boundaries.
      a_neighbors.clear();
      Box neighborsBox(baseCell - IntVect::Unit,
                       baseCell + IntVect::Unit);
      for (BoxIterator bitNbr(neighborsBox); bitNbr.ok(); ++bitNbr)
        {
          IntVect nbr = bitNbr();
          MBStencilElement elt(nbr, baseBlock);
          a_neighbors.push_back(elt);
        }
    }
  else
    { // harder case: baseCell is adjacent to a block boundary.
      Box shiftCellToVertexBox(IntVect::Zero, IntVect::Unit);
      // Fill allVertexElements with all vertex stencil elements that are
      // vertices of cell stencil element a_base.
      Vector<MBStencilElement> allVertexElements;
      for (BoxIterator bitSh(shiftCellToVertexBox); bitSh.ok(); ++bitSh)
        {
          IntVect cellToVertex = bitSh();
          IntVect vertex = baseCell + cellToVertex;
          MBStencilElement elt(vertex, baseBlock);
          // Note that as a NODE, vertex is still valid within baseBlock.
          // First, find all aliases of vertex.
          Vector<MBStencilElement> aliases;
          vertexStencilElementAliases(aliases, elt);
          allVertexElements.append(aliases);
        }
      validStencilElementsFromVertices(a_neighbors, allVertexElements);
    }
  return;
}

void
MultiBlockUtil::vertexAliases(Vector<IntVect>& a_aliasVertices,
                              Vector<int>& a_aliasBlocks,
                              const IntVect& a_baseVertex,
                              int a_baseBlock) const
{
  // It might be better to have a struct consisting of IntVect and int.
  List<IntVect> aliasVerticesList;
  List<int> aliasBlocksList;
  List<IntVect> toVisitVertices;
  List<int> toVisitBlocks;
  toVisitVertices.add(a_baseVertex);
  toVisitBlocks.add(a_baseBlock);
  while (toVisitVertices.isNotEmpty())
    {
      // pop this vertex from toVisit list
      IntVect thisVertex = toVisitVertices.firstElement();
      int thisBlock = toVisitBlocks.firstElement();
      toVisitVertices.removeFirst();
      toVisitBlocks.removeFirst();

      // add this vertex to aliases
      aliasVerticesList.add(thisVertex);
      aliasBlocksList.add(thisBlock);

      Box thisBlockNodesBox = surroundingNodes(m_mappingBlocks[thisBlock]);
      int faceID = 0;
      for (SideIterator sit; sit.ok(); sit.next())
        {
          Side::LoHiSide side = sit();
          IntVect sideBlockEnd = thisBlockNodesBox.sideEnd(side);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (thisVertex[idir] == sideBlockEnd[idir])
                {
                  const BlockBoundary& bb = m_boundaries[thisBlock][faceID];
                  if (bb.isConformal() || bb.isMapped())
                    {
                      IndicesTransformation tfm = bb.getTransformation();
                      int nbrBlock = bb.neighbor();
                      IntVect nbrVertex = tfm.transformNode(thisVertex);
                      if ( !aliasVerticesList.includes(nbrVertex) )
                        if ( !toVisitVertices.includes(nbrVertex) )
                          {
                            toVisitVertices.add(nbrVertex);
                            toVisitBlocks.add(nbrBlock);
                          }
                    }
                }
              faceID++;
            }
        }
    }

  a_aliasVertices.clear();
  a_aliasBlocks.clear();
  while (aliasVerticesList.isNotEmpty())
    {
      a_aliasVertices.push_back(aliasVerticesList.firstElement());
      a_aliasBlocks.push_back(aliasBlocksList.firstElement());
      aliasVerticesList.removeFirst();
      aliasBlocksList.removeFirst();
    }
  return;
}

void
MultiBlockUtil::vertexStencilElementAliases(Vector<MBStencilElement>& a_aliases,
                                            const MBStencilElement& a_base) const
{
  List<MBStencilElement> aliasList;
  List<MBStencilElement> toVisit;
  toVisit.add(a_base);
  while (toVisit.isNotEmpty())
    {
      // pop this stencil element from toVisit list
      const MBStencilElement& thisElement = toVisit.firstElement();
      toVisit.removeFirst();

      // add this stencil element to aliases
      aliasList.add(thisElement);

      int thisBlock = thisElement.block();
      IntVect thisVertex = thisElement.cell();

      Box thisBlockNodesBox = surroundingNodes(m_mappingBlocks[thisBlock]);
      int faceID = 0;
      for (SideIterator sit; sit.ok(); sit.next())
        {
          Side::LoHiSide side = sit();
          IntVect sideBlockEnd = thisBlockNodesBox.sideEnd(side);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (thisVertex[idir] == sideBlockEnd[idir])
                {
                  const BlockBoundary& bb = m_boundaries[thisBlock][faceID];
                  if (bb.isConformal() || bb.isMapped())
                    {
                      IndicesTransformation tfm = bb.getTransformation();
                      int nbrBlock = bb.neighbor();
                      IntVect nbrVertex = tfm.transformNode(thisVertex);
                      MBStencilElement nbr(nbrVertex, nbrBlock);
                      if ( !aliasList.includes(nbr) )
                        if ( !toVisit.includes(nbr) )
                          {
                            toVisit.add(nbr);
                          }
                    }
                }
              faceID++;
            }
        }
    }

  a_aliases.clear();
  while (aliasList.isNotEmpty())
    {
      a_aliases.push_back(aliasList.firstElement());
      aliasList.removeFirst();
    }
  return;
}

void
MultiBlockUtil::validCellsFromVertices(Vector<IntVect>& a_cellIndices,
                                       Vector<int>& a_cellBlocks,
                                       const Vector<IntVect>& a_vertexIndices,
                                       const Vector<int>& a_vertexBlocks) const
{
  Box nodeToCellBox(-IntVect::Unit, IntVect::Zero);
  List<IntVect> cellIndicesList;
  List<int> cellBlocksList;
  for (int i = 0; i < a_vertexIndices.size(); i++)
    {
      IntVect thisVertex = a_vertexIndices[i];
      int thisBlock = a_vertexBlocks[i];
      for (BoxIterator bitOff(nodeToCellBox); bitOff.ok(); ++bitOff)
        {
          IntVect nodeToCell = bitOff();
          IntVect cell = thisVertex + nodeToCell;
          if (m_mappingBlocks[thisBlock].contains(cell))
            if ( !cellIndicesList.includes(cell) )
              {
                cellIndicesList.add(cell);
                cellBlocksList.add(thisBlock);
              }
        }
    }

  a_cellIndices.clear();
  a_cellBlocks.clear();
  while (cellIndicesList.isNotEmpty())
    {
      a_cellIndices.push_back(cellIndicesList.firstElement());
      a_cellBlocks.push_back(cellBlocksList.firstElement());
      cellIndicesList.removeFirst();
      cellBlocksList.removeFirst();
    }
  return;
}

void
MultiBlockUtil::validStencilElementsFromVertices(Vector<MBStencilElement>& a_stencilElements,
                                                 const Vector<MBStencilElement>& a_vertices) const
{
  Box nodeToCellBox(-IntVect::Unit, IntVect::Zero);
  List<MBStencilElement> cellStencilsList;
  for (int i = 0; i < a_vertices.size(); i++)
    {
      const MBStencilElement& thisVertexElement = a_vertices[i];
      const IntVect& thisVertexIndices = thisVertexElement.cell();
      int thisBlock = thisVertexElement.block();
      for (BoxIterator bitOff(nodeToCellBox); bitOff.ok(); ++bitOff)
        {
          IntVect nodeToCell = bitOff();
          IntVect cell = thisVertexIndices + nodeToCell;
          if (m_mappingBlocks[thisBlock].contains(cell))
            { // Add (cell, thisBlock) to list, if it not already in it.
              MBStencilElement cellElement(cell, thisBlock);
              if ( !cellStencilsList.includes(cellElement) )
                {
                  cellStencilsList.add(cellElement);
                }
            }
        }
    }

  a_stencilElements.clear();
  while (cellStencilsList.isNotEmpty())
    {
      a_stencilElements.push_back(cellStencilsList.firstElement());
      cellStencilsList.removeFirst();
    }
  return;
}

bool
MultiBlockUtil::validCellShift(IntVect&        a_shiftedCell,
                               int&            a_shiftedBlock,
                               const IntVect&  a_origCell,
                               int             a_origBlock,
                               int             a_dir,
                               int             a_shift) const
{
  CH_assert(m_isDefined);
  a_shiftedCell = a_origCell;
  Box origBox = m_mappingBlocks[a_origBlock];
  CH_assert(origBox.contains(a_origCell));
  a_shiftedCell.shift(a_dir, a_shift);
  if (origBox.contains(a_shiftedCell))
    { // still in the same block
      a_shiftedBlock = a_origBlock;
      return true;
    }
  else
    {
      int faceID = a_dir;
      if (a_shift > 0) faceID += SpaceDim;

      const BlockBoundary& bb = m_boundaries[a_origBlock][faceID];
      if (bb.isDomainBoundary())
        {
          return false;
        }
      else
        {
          a_shiftedBlock = bb.neighbor();
          IndicesTransformation tfm = bb.getTransformation();
          a_shiftedCell = tfm.transform(a_shiftedCell);
          return true;
        }
    }
}

bool
MultiBlockUtil::validCellShift(MBStencilElement& a_shiftedElt,
                               const MBStencilElement& a_origElt,
                               int             a_dir,
                               int             a_shift) const
{
  CH_assert(m_isDefined);
  const IntVect& origCell = a_origElt.cell();
  int origBlock = a_origElt.block();
  Box origBox = m_mappingBlocks[origBlock];
  CH_assert(origBox.contains(origCell));

  IntVect shiftedCell = origCell;
  shiftedCell.shift(a_dir, a_shift);
  int shiftedBlock = origBlock;
  bool retval = true;
  if (!origBox.contains(shiftedCell))
    {
      int faceID = a_dir;
      if (a_shift > 0) faceID += SpaceDim;

      const BlockBoundary& bb = m_boundaries[origBlock][faceID];
      if (bb.isDomainBoundary())
        {
          retval = false;
        }
      else
        {
          shiftedBlock = bb.neighbor();
          IndicesTransformation tfm = bb.getTransformation();
          shiftedCell = tfm.transform(shiftedCell);
          retval = true;
        }
    }
  a_shiftedElt = MBStencilElement(shiftedCell, shiftedBlock);
  return retval;
}

void
MultiBlockUtil::getValid(IVSFAB<IntVect>& a_validIndices,
                         IVSFAB<int>& a_validBlock,
                         IVSFAB<RealVect>& a_validMappedCenter,
                         const IntVectSet& a_ghostCellsIVS,
                         int a_baseBlockNum) const
{
  CH_assert(m_isDefined);
  const NewCoordSys* baseCoordSysPtr =
    m_coordSysPtr->getCoordSys(a_baseBlockNum);
  for (IVSIterator ivsit(a_ghostCellsIVS); ivsit.ok(); ++ivsit)
    {
      IntVect thisGhostCell = ivsit();
      // baseCenter:  center of thisGhostCell in
      // mapped coordinates of block a_baseBlockNum
      RealVect baseCenter =
        baseCoordSysPtr->centerMappedCoordinates(thisGhostCell);

      // Point baseCenter in baseBlockNum's mapped coordinates
      // becomes validCenter in validBlockNum's mapped coordinates,
      // where validBlockNum is the index of the block where the
      // point is valid.
      int validBlockNum;
      RealVect validCenter;
      m_coordSysPtr->blockRemapping(validCenter, validBlockNum,
                                    baseCenter, a_baseBlockNum);
      CH_assert(validBlockNum >= 0);
      // validCell is in index space of validBlockNum,
      // transformed from
      // thisGhostCell in index space of baseBlockNum.
      IntVect validCell;
      const NewCoordSys* validCoordSysPtr =
        m_coordSysPtr->getCoordSys(validBlockNum);
      const RealVect& validDxi = validCoordSysPtr->dx();
      // WAS for (int idir = 0; idir < SpaceDim; idir++)
      for (int ind = 0; ind < m_interpDimsVect.size(); ind++)
        {
          int idir = m_interpDimsVect[ind];
          validCell[idir] = floor(validCenter[idir] / validDxi[idir]);
        }
      for (int ind = 0; ind < m_fixedDimsVect.size(); ind++)
        {
          int idir = m_fixedDimsVect[ind];
          validCell[idir] = thisGhostCell[idir];
        }
      // You're in trouble if you go outside the domain.
      // So find the closest valid cell that is in the domain.
      const Box& validBox = m_mappingBlocks[validBlockNum];
      if (!validBox.contains(validCell))
        {
          // This should happen only when right on a block boundary.
          // WAS for (int idir = 0; idir < SpaceDim; idir++)
          for (int ind = 0; ind < m_interpDimsVect.size(); ind++)
            {
              int idir = m_interpDimsVect[ind];
              if (validCell[idir] < validBox.smallEnd(idir))
                {
                  validCell[idir] = validBox.smallEnd(idir);
                }
              else if (validCell[idir] > validBox.bigEnd(idir))
                {
                  validCell[idir] = validBox.bigEnd(idir);
                }
            }
        }
      a_validIndices(thisGhostCell, 0) = validCell;
      a_validBlock(thisGhostCell, 0) = validBlockNum;
      a_validMappedCenter(thisGhostCell, 0) = validCenter;
    }
}

void
MultiBlockUtil::getValid(BaseFab<IntVect>& a_validIndices,
                         BaseFab<int>& a_validBlock,
                         BaseFab<RealVect>& a_validMappedCenter,
                         const Box& a_bx,
                         int a_baseBlockNum) const
{
  CH_assert(m_isDefined);
  const NewCoordSys* baseCoordSysPtr =
    m_coordSysPtr->getCoordSys(a_baseBlockNum);
  for (BoxIterator bit(a_bx); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      // baseCenter:  center of iv in
      // mapped coordinates of block a_baseBlockNum
      RealVect baseCenter =
        baseCoordSysPtr->centerMappedCoordinates(iv);

      // Point baseCenter in baseBlockNum's mapped coordinates
      // becomes validCenter in validBlockNum's mapped coordinates,
      // where validBlockNum is the index of the block where the
      // point is valid.
      int validBlockNum;
      RealVect validCenter;
      m_coordSysPtr->blockRemapping(validCenter, validBlockNum,
                                    baseCenter, a_baseBlockNum);
      CH_assert(validBlockNum >= 0);
      // validCell is in index space of validBlockNum,
      // transformed from
      // iv in index space of baseBlockNum.
      IntVect validCell;
      const NewCoordSys* validCoordSysPtr =
        m_coordSysPtr->getCoordSys(validBlockNum);
      const RealVect& validDxi = validCoordSysPtr->dx();
      // WAS for (int idir = 0; idir < SpaceDim; idir++)
      for (int ind = 0; ind < m_interpDimsVect.size(); ind++)
        {
          int idir = m_interpDimsVect[ind];
          validCell[idir] = floor(validCenter[idir] / validDxi[idir]);
        }
      for (int ind = 0; ind < m_fixedDimsVect.size(); ind++)
        {
          int idir = m_fixedDimsVect[ind];
          validCell[idir] = iv[idir];
        }
      // You're in trouble if you go outside the domain.
      // So find the closest valid cell that is in the domain.
      const Box& validBox = m_mappingBlocks[validBlockNum];
      if (!validBox.contains(validCell))
        {
          // This should happen only when right on a block boundary.
          // WAS for (int idir = 0; idir < SpaceDim; idir++)
          for (int ind = 0; ind < m_interpDimsVect.size(); ind++)
            {
              int idir = m_interpDimsVect[ind];
              if (validCell[idir] < validBox.smallEnd(idir))
                {
                  validCell[idir] = validBox.smallEnd(idir);
                }
              else if (validCell[idir] > validBox.bigEnd(idir))
                {
                  validCell[idir] = validBox.bigEnd(idir);
                }
            }
        }
      a_validIndices(iv, 0) = validCell;
      a_validBlock(iv, 0) = validBlockNum;
      a_validMappedCenter(iv, 0) = validCenter;
    }
}

void
MultiBlockUtil::getStencilCells(IVSFAB<MBStencil>& a_stencils,
                                Vector<Box>& a_minBox,
                                const IntVectSet& a_ghostCellsIVS,
                                const IVSFAB<IntVect>& a_validIndices,
                                const IVSFAB<int>& a_validBlock,
                                const IVSFAB<RealVect>& a_validMappedCenter,
                                int a_order,
                                int a_radius) const
{
  CH_assert(m_isDefined);
  // Start with empty a_minBox.
  a_minBox.resize(m_nblocks);
  for (IVSIterator ivsit(a_ghostCellsIVS); ivsit.ok(); ++ivsit)
    {
      IntVect thisGhostCell = ivsit();

      MBStencil stencil;
      // cell in block validBlockNum where validMappedCenterHere is
      const IntVect& validIndices = a_validIndices(thisGhostCell, 0);
      // which valid block contains thisGhostCell
      int validBlockNum = a_validBlock(thisGhostCell, 0);
      // center of thisGhostCell in valid block's mapped coordinates
      const RealVect& validMappedCenterHere = a_validMappedCenter(thisGhostCell, 0);
      getStencilCells(stencil, a_minBox, thisGhostCell,
                      validIndices, validBlockNum, validMappedCenterHere,
                      a_order, a_radius);
      a_stencils(thisGhostCell, 0) = stencil;
    }
}

void
MultiBlockUtil::getStencilCells(MBStencil& a_stencil,
                                Vector<Box>& a_minBox,
                                const IntVect& a_ghostIndices,
                                const IntVect& a_validIndices,
                                int a_validBlockNum,
                                const RealVect& a_validMappedCenter,
                                int a_order,
                                int a_radius) const
{
  MBStencilElement validStencilElement(a_validIndices, a_validBlockNum);

  Vector<MBStencilElement> stencilVec;
  
  if (a_order <= 1)
    {
      stencilVec.push_back(validStencilElement);
    }
  else
    {
      // Set centralCell to be the central cell of the neighborhood,
      // in a_validBlockNum's index space.  It is the cell closest to
      // thisValidCell that is separated from the external boundary
      // by at least a_radius cells.
      // (So centralCell == a_validIndices except near external boundary.)
      // We use only 4th or 5th order, and for those, a_radius == 1.
      IntVect centralCell =
        m_coordSysPtr->cellAvoidingBoundaries(a_validIndices,
                                              a_validBlockNum,
                                              a_radius);
      if (m_fixedDimsVect.size() > 0)
        {
          for (int ind = 0; ind < m_fixedDimsVect.size(); ind++)
            {
              int idir = m_fixedDimsVect[ind];
              centralCell[idir] = a_ghostIndices[idir];
            }
        }

      MBStencilElement baseStencilElement(centralCell, a_validBlockNum);
      // Stencil's inner set:
      // All valid cells having a vertex in common with centralCell.
      commonVertexStencilElements(stencilVec, baseStencilElement);
      if (m_fixedDimsVect.size() > 0)
        { // Remove all stencil cells that do not have 
          // the fixed dimensions set to what they should be.
          Vector<MBStencilElement> stencilVecUpdated;
          int stencilSize = stencilVec.size();
          for (int icell = 0; icell < stencilSize; icell++)
            {
              const MBStencilElement& elt = stencilVec[icell];
              IntVect thisCell = elt.cell();
              bool thisIsIn = true;
              for (int ind = 0; ind < m_fixedDimsVect.size(); ind++)
                {
                  int idir = m_fixedDimsVect[ind];
                  if (thisCell[idir] != a_ghostIndices[idir])
                    {
                      thisIsIn = false;
                    }
                }
              if (thisIsIn)
                {
                  stencilVecUpdated.push_back(elt);
                }
            }
          stencilVec = Vector<MBStencilElement>(stencilVecUpdated);
        }

      if (a_order > 2)
        {
          // To stencilVec, append
          // a_validIndices +/- 2*e, whenever within physical boundary.
          // Then we'll have either 4 or 5 in a row in every dimension.
          // Note this is from a_validIndices, not centralCell.
          IntVect offCenter = centralCell - a_validIndices;
          { CH_TIME("add to stencil");
            for (int ind = 0; ind < m_interpDimsVect.size(); ind++)
              {
                int idir = m_interpDimsVect[ind];
                int offCenterDir = offCenter[idir];
                for (SideIterator sit; sit.ok(); sit.next())
                  {
                    Side::LoHiSide side = sit();
                    int sgn = sign(side);
                    // When offCenterDir == 0, can shift from a_validIndices
                    // by either -2 or +2.
                    
                    // When offCenterDir == -1, we have shifted left
                    // from a_validIndices to centralCell,
                    // so no more room on the right; do iff sgn == -1,
                    // and shift from a_validIndices by -3.
                    
                    // When offCenterDir == +1, we have shifted right
                    // from a_validIndices to centralCell,
                    // so no more room on the left; do iff sgn == +1,
                    // and shift from a_validIndices by +3.
                    if (sgn * offCenterDir >= 0)
                      {
                        int shiftAmount = sgn*2 + offCenterDir;
                        MBStencilElement shiftedElt;
                        if (validCellShift(shiftedElt, validStencilElement,
                                           idir, shiftAmount))
                          {
                            stencilVec.push_back(shiftedElt);
                          }
                      }
                  }
              }
          }
        }
    }

  // So much for a_order == 4.
  // If a_order == 5, need more stencil elements in stencilVec.

  a_stencil = RefCountedPtr< Vector<MBStencilElement> >
    (new Vector<MBStencilElement>(stencilVec));

  MBStencilIterator stencilit(a_stencil);
  for (stencilit.begin(); stencilit.ok(); ++stencilit)
    {
      const MBStencilElement& elt = stencilit();
      IntVect thisStencilCell = elt.cell();
      int thisStencilBlock = elt.block();
      //      CH_assert(m_mappingBlocks[thisStencilBlock].contains(thisStencilCell));
      // Update a_minBox[thisStencilBlock].
      Box& minBoxSrc = a_minBox[thisStencilBlock];
      if ( !minBoxSrc.contains(thisStencilCell) )
        {
          Box thisStencilCellBox(thisStencilCell, thisStencilCell);
          minBoxSrc.minBox(thisStencilCellBox);
        }
    }

  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused++;
}


void
MultiBlockUtil::getVectorStencilCells(IVSFAB<MBVectorStencil>& a_vectorstencils,
                                      const IVSFAB<MBStencil>& a_stencils,
                                      const IntVectSet& a_ghostCellsIVS) const
{
  for (IVSIterator ivsit(a_ghostCellsIVS); ivsit.ok(); ++ivsit)
    {
      IntVect thisGhostCell = ivsit();

      const MBStencil& stencil = a_stencils(thisGhostCell, 0);
      int npts = stencil.size();

      Vector<MBVectorStencilElement> vectorstencilVec(npts);
      for (int ipt = 0; ipt < npts; ipt++)
        {
          const MBStencilElement& stencilElement = stencil[ipt];
          const IntVect& cell = stencilElement.cell();
          int block = stencilElement.block();
          vectorstencilVec[ipt] = MBVectorStencilElement(cell, block);
        }

      MBVectorStencil vectorstencil =
        RefCountedPtr< Vector<MBVectorStencilElement> >
        (new Vector<MBVectorStencilElement>(vectorstencilVec));

      a_vectorstencils(thisGhostCell, 0) = vectorstencil;
    }
}

void
MultiBlockUtil::displacementPowers(Vector<Real>&  a_powers,
                                   const IntVect&     a_cell,
                                   const RealVect&    a_commonMappedBasePoint,
                                   Real               a_avgDistance,
                                   int                a_thisBlockNum,
                                   int                a_commonBlockNum,
                                   const Vector<IntVect>&   a_exponents) const
{ CH_TIME("displacementPowers");
  // This function will return
  // a_powers[pvec] =
  // <((Xvec() - Xvec(a_commonMappedBasePoint))/a_avgDistance)^pvec>_(a_cell),
  // where a_cell is in index space of block a_thisBlockNum.
  const NewCoordSys* coordSysThisBlockPtr =
    m_coordSysPtr->getCoordSys(a_thisBlockNum);

  Vector<RealVect> convCellMappedCenters(m_convSize);
  for (int iconv = 0; iconv < m_convSize; iconv++)
    {
      IntVect convCell = a_cell + m_convBaseCells[iconv];
      convCellMappedCenters[iconv] =
        coordSysThisBlockPtr->centerMappedCoordinates(convCell);
    }

  Vector<int> convSrcBlocks(m_convSize, a_thisBlockNum);

  // For each cell i in the convolution stencil,
  // RealVect convDisplacements[i] is displacement
  // from the point a_commonMappedBasePoint in block a_commonBlockNum
  // to the point convCellMappedCenters[i] in block convSrcBlocks[i].
  Vector<RealVect> convDisplacements =
    m_coordSysPtr->displacements(convCellMappedCenters,
                                 convSrcBlocks,
                                 a_commonMappedBasePoint,
                                 a_commonBlockNum);

  a_powers.assign(0.);
  for (int iconv = 0; iconv < m_convSize; iconv++)
    { CH_TIME("powers of displacement");
      RealVect convDisplacementsNormalized =
        convDisplacements[iconv] / a_avgDistance;

      // Increment a_powers += m_convWeight[iconv] *
      // ((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec_(convCell),
      // where Xvec() has real coordinates,
      // and convCell is in index space of thisStencilCellBlock.
      addPowersPoint(a_powers, convDisplacementsNormalized,
                     a_exponents, m_convWeight[iconv]);
    }
}
      
void
MultiBlockUtil::displacementPowersTransformed(Vector<Real>&  a_powers,
                                              const BaseFab<VectorTransformation>& a_vectorStencilTransformations,
                                              const IntVect&     a_cell,
                                              const RealVect&    a_commonMappedBasePoint,
                                              Real               a_avgDistance,
                                              int                a_thisBlockNum,
                                              int                a_commonBlockNum,
                                              int                a_ghostBlockNum,
                                              const Vector<IntVect>&   a_exponents) const
{ CH_TIME("displacementPowersTransformed");
  // This function will return
  // a_powers[pvec] =
  // <T(Xvec) * ((Xvec() - Xvec(a_commonMappedBasePoint))/a_avgDistance)^pvec>_(a_cell),
  // where a_cell is in index space of block a_thisBlockNum.
  // We do this for all SpaceDim*SpaceDim components of T(Xvec).
  const NewCoordSys* coordSysThisBlockPtr =
    m_coordSysPtr->getCoordSys(a_thisBlockNum);

  Vector<RealVect> convCellMappedCenters(m_convSize); // 1+2*nInterpDims
  for (int iconv = 0; iconv < m_convSize; iconv++)
    {
      IntVect convCell = a_cell + m_convBaseCells[iconv];
      convCellMappedCenters[iconv] =
        coordSysThisBlockPtr->centerMappedCoordinates(convCell);
    }

  Vector<int> convSrcBlocks(m_convSize, a_thisBlockNum);

  // For each cell i in the convolution stencil,
  // RealVect convDisplacements[i] is displacement
  // from point a_commonMappedBasePoint in block a_commonBlockNum
  // to point convCellMappedCenters[i] in block convSrcBlocks[i].
  Vector<RealVect> convDisplacements =
    m_coordSysPtr->displacements(convCellMappedCenters,
                                 convSrcBlocks,
                                 a_commonMappedBasePoint,
                                 a_commonBlockNum);

  // int numTaylorCoeffs = a_powers.size() / (SpaceDim*SpaceDim);
  a_powers.assign(0.); // length SpaceDim*SpaceDim * numTaylorCoeffs

  for (int iconv = 0; iconv < m_convSize; iconv++)
    { CH_TIME("powers of displacement");
      RealVect convDisplacementsNormalized =
        convDisplacements[iconv] / a_avgDistance;

      // Set basePowers = m_convWeight[iconv] *
      // ((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec_(convCell),
      // where Xvec() has real coordinates,
      // and convCell is in index space of thisStencilCellBlock.
      Real weight = m_convWeight[iconv];
      // Multiply weight by T(Xvec), where we get T(XVec) from
      // m_coordSysPtr function that uses
      // convCellMappedCenters, convSrcBlocks, a_commonBlockNum.

      // For each cell i in the convolution stencil,
      // RealVect convDisplacements[i] is displacement
      // from point a_commonMappedBasePoint in block a_commonBlockNum
      // to point convCellMappedCenters[i] in block convSrcBlocks[i].
      
      // VectorTransformation vectorBlockTransformation(int a_nDst,
      // const RealVect& a_xiSrc, int a_nSrc) const;
      // VectorTransformation vtToGhost transforms components of a vector
      // in basis of block thisVectorStencilCellBlock
      // to components in a basis in block a_ghostBlockNum.
      // The vector is situated at point thisVectorStencilCellCenter
      // of block thisVectorStencilCellBlock.
      IntVect convCell = a_cell + m_convBaseCells[iconv];
      VectorTransformation vt = a_vectorStencilTransformations(convCell, 0);

      // Now set a_powers += basePowers * T(Xvec).
      addPowersVector(a_powers, convDisplacementsNormalized,
                      a_exponents, weight, vt.dataPtr());
    }
}
      
void
MultiBlockUtil::getVectorTransformations(Vector< BaseFab<VectorTransformation>* >& a_vectorTransformations,
                                         const IVSFAB<MBVectorStencil>& a_vectorstencilsFab,
                                         const IndexType& a_type,
                                         const IntVectSet& a_ghostCellsIVS,
                                         int a_ghostBlockNum) const
{ CH_TIME("MultiBlockUtil::getVectorTransformations");
  Vector< BaseFab<bool>* > gotVectorTransformations(m_nblocks);
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    {
      const Box& bxSrc = a_vectorTransformations[srcBlock]->box();
      gotVectorTransformations[srcBlock] = new BaseFab<bool>(bxSrc, 1);
      gotVectorTransformations[srcBlock]->setVal(false);
    }
  for (IVSIterator ivsit(a_ghostCellsIVS); ivsit.ok(); ++ivsit)
    {
      IntVect thisGhostCell = ivsit();
      const MBVectorStencil& thisGhostStencil = a_vectorstencilsFab(thisGhostCell, 0);
      int vectorstencilSize = thisGhostStencil.size();
      for (int icell = 0; icell < vectorstencilSize; icell++)
        {
          const MBVectorStencilElement& vectorstencilElement =
            thisGhostStencil[icell];
          const IntVect& thisVectorStencilCellIndices =
            vectorstencilElement.cell();
          int thisVectorStencilCellBlock =
            vectorstencilElement.block();
          BaseFab<VectorTransformation>& vectorTransformationsSrc =
            *(a_vectorTransformations[thisVectorStencilCellBlock]);
          BaseFab<bool>& gotVectorTransformationsSrc =
            *(gotVectorTransformations[thisVectorStencilCellBlock]);
          if (a_type == IndexType::TheNodeType())
            {
              if ( !gotVectorTransformationsSrc(thisVectorStencilCellIndices, 0) )
                {
                  VectorTransformation vtToGhost =
                    m_coordSysPtr->vectorBlockTransformationCenter(a_ghostBlockNum,
                                                                   thisVectorStencilCellIndices,
                                                                   thisVectorStencilCellBlock);
                  
                  vectorTransformationsSrc(thisVectorStencilCellIndices, 0) =
                    vtToGhost.inverse();
                  gotVectorTransformationsSrc(thisVectorStencilCellIndices, 0) = true;
                }
            }
          else if (a_type == IndexType::TheCellType())
            {
              for (int iconv = 0; iconv < m_convSize; iconv++)
                {
                  IntVect convCellIndices = 
                    thisVectorStencilCellIndices + m_convBaseCells[iconv];
                  if ( !gotVectorTransformationsSrc(convCellIndices, 0) )
                    {
                      VectorTransformation vtToGhost =
                        m_coordSysPtr->vectorBlockTransformationCenter(a_ghostBlockNum,
                                                                       convCellIndices,
                                                                       thisVectorStencilCellBlock);
                      vectorTransformationsSrc(convCellIndices, 0) =
                        vtToGhost.inverse();
                      gotVectorTransformationsSrc(convCellIndices, 0) = true;
                    }
                }
            }
        }
    }
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    {
      delete gotVectorTransformations[srcBlock];
    }
}


void
MultiBlockUtil::getWeights(IVSFAB<MBStencil>& a_stencilsFab,
                           const IndexType& a_type,
                           const IntVectSet& a_ghostCellsIVS,
                           const IVSFAB<IntVect>& a_validIndices,
                           const IVSFAB<int>& a_validBlock,
                           const IVSFAB<RealVect>& a_validMappedCenter,
                           const Vector<IntVect>& a_exponents,
                           int a_ghostBlockNum) const
{ CH_TIME("MultiBlockUtil::getWeights");
  CH_assert(m_isDefined);
  for (IVSIterator ivsit(a_ghostCellsIVS); ivsit.ok(); ++ivsit)
    {
      IntVect thisGhostCell = ivsit();
      MBStencil& thisGhostStencil = a_stencilsFab(thisGhostCell, 0);
      // cell in block validBlockNum where validMappedCenterHere is
      const IntVect& validIndices = a_validIndices(thisGhostCell, 0);
      // which valid block contains thisGhostCell
      int validBlockNum = a_validBlock(thisGhostCell, 0);
      // center of thisGhostCell in valid block's mapped coordinates
      const RealVect& validMappedCenterHere = a_validMappedCenter(thisGhostCell, 0);
      getWeights(thisGhostStencil, a_type, thisGhostCell,
                 validIndices, validBlockNum, validMappedCenterHere,
                 a_exponents, a_ghostBlockNum);
    }
}

void
MultiBlockUtil::getWeights(MBStencil& a_stencil,
                           const IndexType& a_type,
                           const IntVect& a_ghostIndices,
                           const IntVect& a_validIndices,
                           int a_validBlockNum,
                           const RealVect& a_validMappedCenter, 
                           const Vector<IntVect>& a_exponents,
                           int a_ghostBlockNum) const
{ CH_TIME("MultiBlockUtil::getWeights on one");
  int stencilSize = a_stencil.size();

  // Set to total of the distances, in real space,
  // from center of ghost cell to center of each stencil cell.
  Real totalDistance = 0.;
  // Mapped-space coordinates of center of each valid cell in the stencil.
  Vector<RealVect> stencilCellCenters(stencilSize);
  // Block of each valid cell in the stencil.
  Vector<int> stencilBlocks(stencilSize);
  { CH_TIME("total distance");
    for (int icell = 0; icell < stencilSize; icell++)
      {
        const MBStencilElement& stencilElement = a_stencil[icell];
        const IntVect& thisStencilCellIndices = stencilElement.cell();
        int thisStencilCellBlock = stencilElement.block();
        stencilBlocks[icell] = thisStencilCellBlock;
        stencilCellCenters[icell] =
          m_coordSysPtr->getCoordSys(thisStencilCellBlock)
          ->centerMappedCoordinates(thisStencilCellIndices);
        // Center of ghost cell in mapped coordinates is
        // a_validMappedCenter in block a_validBlockNum.
        totalDistance +=
          m_coordSysPtr->distance(stencilCellCenters[icell],
                                  thisStencilCellBlock,
                                  a_validMappedCenter,
                                  a_validBlockNum);
      }
  }
  Real avgDistance = totalDistance / Real(stencilSize);

  // a_exponents contains the exponents, SpaceDim at a time.
  // With 2D and degree 3 and no fixed dimensions, these are:
  // 00 10 20 30 01 11 21 02 12 03.
  // These are Vector<int> instead of Vector<IntVect> because
  // we need to send them to Chombo Fortran.
  int numTaylorCoeffs = a_exponents.size();

  /*
    This is used only if cell-centered.
  */
  Vector<RealVect> disp;
  if (a_type == IndexType::TheNodeType())
    { // cell-centered only
      // Return Vector<RealVect> disp, where for each stencil cell i,
      // RealVect disp[i] is displacement
      // from the point a_validMappedCenter in block a_validBlockNum 
      // to the point vectorstencilCellCenters[i] in block stencilBlocks[i].
      // These are displacements in real space
      // from the center of the ghost cell to the center of each stencil cell.
      // Center of ghost cell in mapped coordinates is
      // a_validMappedCenter in block a_validBlockNum.
      disp = m_coordSysPtr->displacements(stencilCellCenters,
                                          stencilBlocks,
                                          a_validMappedCenter,
                                          a_validBlockNum);
    }

  /*
    If cell-averaged:

    thisDispPowers[pvec] =
    <((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec>_g
    where Xvec() has real coordinates.

    If we're interpolating cell centers instead of cell averages,
    then we don't need thisDispPowers, because
    the only Xvec() will be Xvec(ghostcenter), and
    ((Xvec(ghostcenter) - Xvec(ghostcenter))/avgDistance)^pvec_g
    == delta(pvec, 0).
    In other words, when interpolating cell centers, we use
    only the first coefficient of the Taylor series, and
    ignore the other coefficients.
  */
  Vector<Real> thisDispPowers;
  if (a_type == IndexType::TheCellType())
    { CH_TIME("center displacement powers"); // cell-averaged only
      thisDispPowers.resize(numTaylorCoeffs);
      displacementPowers(thisDispPowers,
                         a_ghostIndices,
                         a_validMappedCenter,
                         avgDistance,
                         a_ghostBlockNum,
                         a_validBlockNum,
                         a_exponents);
    }

  // To hold displacement powers (centered or averaged) of stencil cells,
  // numTaylorCoeffs at a time.
  Vector<Real> coeffs;
  if (a_type == IndexType::TheNodeType())
    { CH_TIME("centered displacement powers"); // cell-centered only
      for (int icell = 0; icell < stencilSize; icell++)
        {
          RealVect dispSrcCell = disp[icell] / avgDistance;
          // Set coeffsCell[pvec] =
          // ((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec_(thisStencilCellIndices),
          // where Xvec() has real coordinates,
          // a_exponents holds the powers pvec,
          // and thisStencilCellIndices in index space of thisStencilCellBlock.
          Vector<Real> coeffsCell(numTaylorCoeffs, 0.);
          addPowersPoint(coeffsCell, dispSrcCell,
                         a_exponents, 1.);
          coeffs.append(coeffsCell);
        } // end loop over cells in stencil
    }
  else if (a_type == IndexType::TheCellType())
    { CH_TIME("averaged displacement powers"); // cell-averaged only
      for (int icell = 0; icell < stencilSize; icell++)
        {
          const MBStencilElement& stencilElement = a_stencil[icell];
          const IntVect& thisStencilCellIndices = stencilElement.cell();
          int thisStencilCellBlock = stencilElement.block();
          
          Vector<Real> dispAvgPowers(numTaylorCoeffs, 0.);
          displacementPowers(dispAvgPowers,
                             thisStencilCellIndices,
                             a_validMappedCenter,
                             avgDistance,
                             thisStencilCellBlock,
                             a_validBlockNum,
                             a_exponents);
          coeffs.append(dispAvgPowers);
        } // end loop over cells in stencil
    }
  else
    {
      MayDay::Error("Bad index type");
    }

  // Now using these powers, find weights.

  // stencilSize equations or function evaluations.
  // numTaylorCoeffs variables as displacement powers or coefficients.

  // Vector<Real> coeffs has length stencilSize * numTaylorCoeffs.
  // - On entry, coeffs holds powers of displacement,
  // numTaylorCoeffs at a time for each stencil cell.
  // - On exit, coeffs is replaced by coefficients,
  // numTaylorCoeffs at a time for each stencil cell.
  LAPACKMatrix A(stencilSize, numTaylorCoeffs);
  {
    int icoeff = 0;
    for (int icell = 0; icell < stencilSize; icell++)
      {
        for (int ico = 0; ico < numTaylorCoeffs; ico++)
          {
            A(icell, ico) = coeffs[icoeff];
            icoeff++;
          }
      }
  }
  A.pseudoInvertUsingQR();
  {
    int icoeff = 0;
    for (int icell = 0; icell < stencilSize; icell++)
      {
        for (int ico = 0; ico < numTaylorCoeffs; ico++)
          {
            coeffs[icoeff] = A(ico, icell);
            icoeff++;
          }
      }
  }

  int stencilBaseIndex = 0;
  for (int icell = 0; icell < stencilSize; icell++)
    {
      Real wt;
      // From coeffs(numTaylorCoeffs*icell+(0:numTaylorCoeffs-1))
      // need to find Real wt,
      // and then do a_stencil[icell].setWeight(wt);
      // This is the weight in the interpolation stencil.
      if (a_type == IndexType::TheNodeType())
        { // cell-centered only:
          // At point Xvec in the ghost cell,
          // weight of icell in thisStencilCellIndices is
          // coeffs[icell, pvec] *
          // ((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec_icell.
          // At cell center, Xvec() = Xvec(ghostcenter),
          // so we can ignore the terms with nonzero pvec, and
          // this is simply coeffs[icell, 0].
          wt = coeffs[stencilBaseIndex];
        }
      else if (a_type == IndexType::TheCellType())
        { // cell-averaged only:
          // At point Xvec in the ghost cell,
          // weight of icell in thisStencilCellIndices is
          // coeffs[icell, pvec] *
          // <((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec>_icell
          // where the average is over the ghost cell.
          wt = dotSubvectors(coeffs, stencilBaseIndex,
                             thisDispPowers, 0,
                             numTaylorCoeffs);
        }
      a_stencil[icell].setWeight(wt);
      stencilBaseIndex += numTaylorCoeffs;
    }

  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused++;
}


void
MultiBlockUtil::getVectorWeights(IVSFAB<MBVectorStencil>& a_vectorstencilsFab,
                                 const Vector< BaseFab<VectorTransformation>* >& a_vectorStencilTransformations,
                                 const IndexType& a_type,
                                 const IntVectSet& a_ghostCellsIVS,
                                 const IVSFAB<IntVect>& a_validIndices,
                                 const IVSFAB<int>& a_validBlock,
                                 const IVSFAB<RealVect>& a_validMappedCenter,
                                 const Vector<IntVect>& a_exponents,
                                 int a_ghostBlockNum) const
{ CH_TIME("MultiBlockUtil::getVectorWeights");
  CH_assert(m_isDefined);
  for (IVSIterator ivsit(a_ghostCellsIVS); ivsit.ok(); ++ivsit)
    {
      IntVect thisGhostCell = ivsit();
      MBVectorStencil& thisGhostVectorStencil =
        a_vectorstencilsFab(thisGhostCell, 0);
      // cell in block validBlockNum where validMappedCenterHere is
      const IntVect& validIndices = a_validIndices(thisGhostCell, 0);
      // which valid block contains thisGhostCell
      int validBlockNum = a_validBlock(thisGhostCell, 0);
      // center of thisGhostCell in valid block's mapped coordinates
      const RealVect& validMappedCenterHere = a_validMappedCenter(thisGhostCell, 0);
      getVectorWeights(thisGhostVectorStencil, a_vectorStencilTransformations,
                       a_type, thisGhostCell,
                       validIndices, validBlockNum, validMappedCenterHere,
                       a_exponents, a_ghostBlockNum);
    }
}


void
MultiBlockUtil::getVectorWeights(MBVectorStencil& a_vectorstencil,
                                 const Vector< BaseFab<VectorTransformation>* >& a_vectorStencilTransformations,
                                 const IndexType& a_type,
                                 const IntVect& a_ghostIndices,
                                 const IntVect& a_validIndices,
                                 int a_validBlockNum,
                                 const RealVect& a_validMappedCenter,
                                 const Vector<IntVect>& a_exponents,
                                 int a_ghostBlockNum) const
{ CH_TIME("MultiBlockUtil::getVectorWeights on one");
  /*
    FIXME: this is a copy of a section of getWeights().
   */
  int vectorstencilSize = a_vectorstencil.size();

  // Set to total of the distances, in real space,
  // from center of ghost cell to center of each stencil cell.
  Real totalDistance = 0.;
  // Indices of each valid cell in the stencil.
  Vector<IntVect> vectorstencilCellIndices(vectorstencilSize);
  // Mapped-space coordinates of center of each valid cell in the stencil.
  Vector<RealVect> vectorstencilCellCenters(vectorstencilSize);
  // Block of each valid cell in the stencil.
  Vector<int> vectorstencilBlocks(vectorstencilSize);
  { CH_TIME("total distance");
    for (int icell = 0; icell < vectorstencilSize; icell++)
      {
        const MBVectorStencilElement& vectorstencilElement = a_vectorstencil[icell];
        const IntVect& thisVectorStencilCellIndices = vectorstencilElement.cell();
        vectorstencilCellIndices[icell] = thisVectorStencilCellIndices;
        int thisVectorStencilCellBlock = vectorstencilElement.block();
        vectorstencilBlocks[icell] = thisVectorStencilCellBlock;
        vectorstencilCellCenters[icell] =
          m_coordSysPtr->getCoordSys(thisVectorStencilCellBlock)
          ->centerMappedCoordinates(thisVectorStencilCellIndices);
        totalDistance +=
          m_coordSysPtr->distance(vectorstencilCellCenters[icell],
                                  thisVectorStencilCellBlock,
                                  a_validMappedCenter,
                                  a_validBlockNum);
      }
  }
  Real avgDistance = totalDistance / Real(vectorstencilSize);

  // a_exponents contains the exponents, SpaceDim at a time.
  // With 2D and degree 3 and no fixed dimensions, these are:
  // 00 10 20 30 01 11 21 02 12 03.
  // These are Vector<int> instead of Vector<IntVect> because
  // we need to send them to Chombo Fortran.
  int numTaylorCoeffs = a_exponents.size();
  int numVectorAllCoeffs = SpaceDim*SpaceDim * numTaylorCoeffs;

  /*
    This is used only if cell-centered.
    FIXME: this is a copy of a section of getWeights().
  */
  Vector<RealVect> disp;
  if (a_type == IndexType::TheNodeType())
    { // cell-centered only
      // Return Vector<RealVect> disp, where for each stencil cell i,
      // RealVect disp[i] is displacement
      // from the point a_validMappedCenter in block a_validBlockNum 
      // to the point vectorstencilCellCenters[i] in block vectorstencilBlocks[i].
      // These are displacements in real space
      // from the center of the ghost cell to the center of each stencil cell.
      // Center of ghost cell in mapped coordinates is
      // a_validMappedCenter in block a_validBlockNum.
      disp = m_coordSysPtr->displacements(vectorstencilCellCenters,
                                          vectorstencilBlocks,
                                          a_validMappedCenter,
                                          a_validBlockNum);
      // For each i, we will get monomials from disp[i],
      // and multiply the monomials by all the components of the
      // basis transform matrix.
    }

  /*
    If cell-averaged:

    thisDispPowersTransformed[dst, src, pvec] =
    T(Xvec())[dst,src] * <((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec>_g
    where Xvec() has real coordinates.

    If we're interpolating cell centers instead of cell averages,
    then we don't need thisDispPowersTransforemd, because
    the only Xvec() will be Xvec(ghostcenter), and
    ((Xvec(ghostcenter) - Xvec(ghostcenter))/avgDistance)^pvec_g
    == delta(pvec, 0).
    In other words, when interpolating cell centers, we use
    only the first coefficient of the Taylor series, and
    ignore the other coefficients.
  */
  Vector<Real> thisDispPowers;
  if (a_type == IndexType::TheCellType())
    { CH_TIME("center displacement powers"); // cell-averaged only
      thisDispPowers.resize(numTaylorCoeffs);
      displacementPowers(thisDispPowers,
                         a_ghostIndices,
                         a_validMappedCenter,
                         avgDistance,
                         a_ghostBlockNum,
                         a_validBlockNum,
                         a_exponents);
    }

  // To hold transformed displacement powers (centered or averaged) of stencil cells,
  // numVectorAllCoeffs = SpaceDim*SpaceDim * numTaylorCoeffs at a time.
  // For each stencil cell, have SpaceDim*SpaceDim groups, each group
  // having numTaylorCoeffs at a time corresponding to a matrix element.
  Vector<Real> coeffs;
  if (a_type == IndexType::TheNodeType())
    { CH_TIME("centered displacement powers"); // cell-centered only
      for (int icell = 0; icell < vectorstencilSize; icell++)
        {
          // Recall: disp[icell] is displacement in real space
          // from the center of the valid cell containing the ghost cell
          // to the center of stencil cell icell.
          RealVect dispSrcCell = disp[icell] / avgDistance;
          // Set baseCoeffsCell[pvec] =
          // ((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec_(thisVectorStencilCellIndices),
          // where Xvec() has real coordinates,
          // a_exponents holds the powers pvec,
          // and thisVectorStencilCellIndices in index space of thisVectorStencilCellBlock.
          Vector<Real> baseCoeffsCell(numTaylorCoeffs, 0.);
          addPowersPoint(baseCoeffsCell, dispSrcCell,
                         a_exponents, 1.);

          // Now set coeffs = baseCoeffsCell * T(Xvec).
          // Need vectorstencilCellCenters[icell] and vectorstencilBlocks[icell].

          // Store in this order for each stencil cell (say 2D, degree 3):
          // [1, x, x^2, x^3, y, x*y, x*2*y, y^2, x*y^2, y^3]*T00,
          // [1, x, x^2, x^3, y, x*y, x*2*y, y^2, x*y^2, y^3]*T01,
          // [1, x, x^2, x^3, y, x*y, x*2*y, y^2, x*y^2, y^3]*T10,
          // [1, x, x^2, x^3, y, x*y, x*2*y, y^2, x*y^2, y^3]*T11.
          // The first two rows above make up one row in the matrix:
          // coefficients of 1st function component evaluated at stencil cell.
          // The second two rows above make up one row in the matrix:
          // coefficients of 2nd function component evaluated at stencil cell.
          // Here (x,y) is the normalized displacement in real space
          // from the ghost cell center to the stencil cell center.
          // ( T00 T01 )
          // ( T10 T11 ) is the transformation matrix at stencil cell center
          // from ghost cell block to stencil cell block.
          // The overdetermined system is (P < N):
          // ( x1^p*y1^q*T00(x1,y1) x1^p*y1^q*T01(x1,y1) )( a1 )   ( f(x1,y1) )
          // ( x1^p*y1^q*T10(x1,y1) x1^p*y1^q*T11(x1,y1) )( :  )   ( g(x1,y1) )
          // (          :                    :           )( aP ) = (  :       )
          // (          :                    :           )( b1 ) = (  :       )
          // ( xN^p*yN^q*T00(xN,yN) xN^p*yN^q*T01(xN,yN) )( :  )   ( f(xN,yN) )
          // ( xN^p*yN^q*T10(xN,yN) xN^p*yN^q*T11(xN,yN) )( bP )   ( g(xN,yN) )
          // because for each stencil cell i,
          // f(xi,yi) = sum_{p,q}
          //   a_{p,q}*xi^p*yi^q*T00(xi,yi) + b_{p,q}*xi^p*yi^q*T01(xi,yi) ;
          // g(xi,yi) = sum_{p,q}
          //   a_{p,q}*xi^p*yi^q*T10(xi,yi) + b_{p,q}*xi^p*yi^q*T11(xi,yi) .
          RealVect thisVectorStencilCellCenter = vectorstencilCellCenters[icell];

          int thisVectorStencilCellBlock = vectorstencilBlocks[icell];
          // VectorTransformation vectorBlockTransformation(int a_nDst,
          // const RealVect& a_xiSrc, int a_nSrc) const;
          // VectorTransformation vtToGhost transforms components of a vector
          // in basis of block thisVectorStencilCellBlock
          // to components in a basis in block a_ghostBlockNum.
          // The vector is situated at point thisVectorStencilCellCenter
          // of block thisVectorStencilCellBlock.
          IntVect thisVectorStencilCellIndices = vectorstencilCellIndices[icell];
          VectorTransformation vt =
            (*a_vectorStencilTransformations[thisVectorStencilCellBlock])(thisVectorStencilCellIndices, 0);
          for (int dst = 0; dst < SpaceDim; dst++)
            {
              for (int src = 0; src < SpaceDim; src++)
                {
                  Real dstTransComp = vt.component(dst, src);
                  Vector<Real> transformedCoeffsCell(numTaylorCoeffs);
                  for (int ico = 0; ico < numTaylorCoeffs; ico++)
                    {
                      transformedCoeffsCell[ico] = dstTransComp *
                        baseCoeffsCell[ico];
                    }
                  coeffs.append(transformedCoeffsCell);
                }
            }
        } // end loop over cells in stencil
    }
  else if (a_type == IndexType::TheCellType())
    { CH_TIME("averaged displacement powers"); // cell-averaged only
      for (int icell = 0; icell < vectorstencilSize; icell++)
        {
          const MBVectorStencilElement& vectorstencilElement = a_vectorstencil[icell];
          const IntVect& thisVectorStencilCellIndices = vectorstencilElement.cell();
          int thisVectorStencilCellBlock = vectorstencilElement.block();
          
          Vector<Real> dispAvgPowers(numVectorAllCoeffs, 0.);
          // There are SpaceDim*SpaceDim components of T(x,y).
          displacementPowersTransformed(dispAvgPowers,
                                        *(a_vectorStencilTransformations[thisVectorStencilCellBlock]),
                                        thisVectorStencilCellIndices,
                                        a_validMappedCenter,
                                        avgDistance,
                                        thisVectorStencilCellBlock,
                                        a_validBlockNum,
                                        a_ghostBlockNum,
                                        a_exponents);
          coeffs.append(dispAvgPowers);
        } // end loop over cells in stencil
    }
  else
    {
      MayDay::Error("Bad index type");
    }

  // Now using these powers, find weights.

  // numFunctionValues = SpaceDim * vectorstencilSize equations or
  // function evaluations.
  // numVectorCoeffs = SpaceDim * numTaylorCoeffs variables as
  // displacement powers or coefficients.

  // Vector<Real> coeffs has length vectorstencilSize * numVectorAllCoeffs.
  // - On entry, coeffs holds a row of transformed powers of displacement,
  // numVectorCoeffs = SpaceDim * numTaylorCoeffs at a time
  // for each function evaluation, and within each function evaluation,
  // grouped by matrix column, numTaylorCoeffs at a time.
  // - On exit, coeffs is replaced by coefficients,
  // numVectorCoeffs = SpaceDim * numTaylorCoeffs at a time
  // for each function evaluation, and within each function evaluation,
  // grouped by matrix column, numTaylorCoeffs at a time.
  int numFunctionValues = SpaceDim * vectorstencilSize;
  int numVectorCoeffs = SpaceDim * numTaylorCoeffs;
  LAPACKMatrix A(numFunctionValues, numVectorCoeffs);
  {
    int icoeff = 0;
    for (int icell = 0; icell < numFunctionValues; icell++)
      {
        for (int ico = 0; ico < numVectorCoeffs; ico++)
          {
            A(icell, ico) = coeffs[icoeff];
            icoeff++;
          }
      }
  }
  A.pseudoInvertUsingQR();
  {
    int icoeff = 0;
    for (int icell = 0; icell < numFunctionValues; icell++)
      {
        for (int ico = 0; ico < numVectorCoeffs; ico++)
          {
            coeffs[icoeff] = A(ico, icell);
            icoeff++;
          }
      }
  }

  int vectorstencilBaseIndex = 0;
  for (int icell = 0; icell < vectorstencilSize; icell++)
    {
      Tuple<Real, SpaceDim*SpaceDim> wt;
      // From coeffs(numVectorAllCoeffs*icell+(0:numVectorAllCoeffs-1))
      // need to find Tuple<Real, SpaceDim*SpaceDim> wt,
      // and then do a_vectorstencil[icell].setWeight(wt);
      // This is the weight in the interpolation stencil.
      int tupleIndex = 0;
      for (int dst = 0; dst < SpaceDim; dst++)
        {
          for (int src = 0; src < SpaceDim; src++)
            {
              // tupleIndex = dst*SpaceDim + src
              if (a_type == IndexType::TheNodeType())
                { // cell-centered only:
                  // At point Xvec in the ghost cell,
                  // weight of icell in thisVectorStencilCellIndices is
                  // coeffs[icell, pvec] *
                  // ((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec_icell.
                  // At cell center, Xvec() = Xvec(ghostcenter),
                  // so we can ignore the terms with nonzero pvec, and
                  // this is simply coeffs[icell, 0].
                  // In 2D:
                  // wt[0] is weight of src=0 for dst=0
                  // wt[1] is weight of src=1 for dst=0
                  // wt[2] is weight of src=0 for dst=1
                  // wt[3] is weight of src=1 for dst=1
                  wt[tupleIndex] = coeffs[vectorstencilBaseIndex];
                }
              else if (a_type == IndexType::TheCellType())
                { // cell-averaged only:
                  // At point Xvec in the ghost cell,
                  // weight of icell in thisVectorStencilCellIndices is
                  // coeffs[icell, pvec] *
                  // <((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec>_icell
                  // where the average is over the ghost cell.
                  wt[tupleIndex] =
                    dotSubvectors(coeffs, vectorstencilBaseIndex,
                                  thisDispPowers, 0,
                                  numTaylorCoeffs);
                  // If I set wt[tupleIndex] = coeffs[vectorstencilBaseIndex];
                  // then I get 2nd-order convergence.
                }
              tupleIndex++;
              vectorstencilBaseIndex += numTaylorCoeffs;
            } // end loop over source components
        } // end loop over dest components
      a_vectorstencil[icell].setWeight(wt);
      // numVectorAllCoeffs = SpaceDim*SpaceDim * numTaylorCoeffs
      // vectorstencilBaseIndex += numVectorAllCoeffs;
    } // end loop over all cells in the stencil

  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused++;
}


void
MultiBlockUtil::copyStencilsFromBlock(IVSFAB<MBStencil>& a_stencilsFab,
                                      const IntVectSet& a_ghostCellsIVS,
                                      const IVSFAB<MBStencil>& a_blockStencilsFab) const
{
  CH_assert(m_isDefined);
  for (IVSIterator ivsit(a_ghostCellsIVS); ivsit.ok(); ++ivsit)
    {
      // on the collapsed layout
      IntVect thisGhostCell = ivsit();
      IntVect fixedCell(thisGhostCell);
      for (int ind = 0; ind < m_fixedDimsVect.size(); ind++)
        {
          int idir = m_fixedDimsVect[ind];
          fixedCell[idir] = m_fixedPt[ind];
        }
      copyStencilsFromBlock(a_stencilsFab(thisGhostCell, 0),
                            thisGhostCell,
                            a_blockStencilsFab(fixedCell, 0));
    }
}      


void
MultiBlockUtil::copyStencilsFromBlock(MBStencil& a_stencil,
                                      const IntVect& a_iv,
                                      const MBStencil& a_blockStencil) const
{
  // a_stencil will be the same as a_blockStencil except that
  // cell[m_fixedDimsVect] for each stencil element will be set to
  // a_iv[m_fixedDimsVect] instead of m_fixedPt.
  Vector<MBStencilElement> stencilVec;
  MBStencilIterator stencilit(a_blockStencil);
  for (stencilit.begin(); stencilit.ok(); ++stencilit)
    {
      const MBStencilElement& blockStencilElement = stencilit();
      const IntVect& blockStencilCell = blockStencilElement.cell();
      int block = blockStencilElement.block(); // keep this
      Real weight = blockStencilElement.weight(); // keep this
      IntVect thisStencilCell = blockStencilCell; // modify this in m_fixedDims
      for (int idir = m_fixedDims.begin(); idir <= m_fixedDims.end(); idir++)
        {
          thisStencilCell[idir] = a_iv[idir];
        }
      MBStencilElement elt(thisStencilCell, block, weight);
      stencilVec.push_back(elt);
    }
  a_stencil = RefCountedPtr< Vector<MBStencilElement> >
    (new Vector<MBStencilElement>(stencilVec));
}


void
MultiBlockUtil::copyVectorStencilsFromBlock(IVSFAB<MBVectorStencil>& a_vectorstencilsFab,
                                            const IntVectSet& a_ghostCellsIVS,
                                            const IVSFAB<MBVectorStencil>& a_blockVectorStencilsFab) const
{
  CH_assert(m_isDefined);
  for (IVSIterator ivsit(a_ghostCellsIVS); ivsit.ok(); ++ivsit)
    {
      // on the collapsed layout
      IntVect thisGhostCell = ivsit();
      IntVect fixedCell(thisGhostCell);
      for (int ind = 0; ind < m_fixedDimsVect.size(); ind++)
        {
          int idir = m_fixedDimsVect[ind];
          fixedCell[idir] = m_fixedPt[ind];
        }
      copyVectorStencilsFromBlock(a_vectorstencilsFab(thisGhostCell, 0),
                                  thisGhostCell,
                                  a_blockVectorStencilsFab(fixedCell, 0));
    }
}      


void
MultiBlockUtil::copyVectorStencilsFromBlock(MBVectorStencil& a_vectorstencil,
                                            const IntVect& a_iv,
                                            const MBVectorStencil& a_blockVectorStencil) const
{
  // a_vectorstencil will be the same as a_blockVectorStencil except that
  // cell[m_fixedDimsVect] for each vector stencil element will be set to
  // a_iv[m_fixedDimsVect] instead of m_fixedPt.
  Vector<MBVectorStencilElement> vectorstencilVec;
  MBVectorStencilIterator vectorstencilit(a_blockVectorStencil);
  for (vectorstencilit.begin(); vectorstencilit.ok(); ++vectorstencilit)
    {
      const MBVectorStencilElement& blockVectorStencilElement = vectorstencilit();
      const IntVect& blockVectorStencilCell = blockVectorStencilElement.cell();
      int block = blockVectorStencilElement.block(); // keep this
      Tuple<Real, SpaceDim*SpaceDim> weight = blockVectorStencilElement.weight(); // keep this
      IntVect thisVectorStencilCell = blockVectorStencilCell; // modify this in m_fixedDims
      for (int idir = m_fixedDims.begin(); idir <= m_fixedDims.end(); idir++)
        {
          thisVectorStencilCell[idir] = a_iv[idir];
        }
      MBVectorStencilElement elt(thisVectorStencilCell, block, weight);
      vectorstencilVec.push_back(elt);
    }
  a_vectorstencil = RefCountedPtr< Vector<MBVectorStencilElement> >
    (new Vector<MBVectorStencilElement>(vectorstencilVec));
}


Box
MultiBlockUtil::boxFixed(const Box& a_bx)
{
  Box returnBox(a_bx);
  for (int ind = 0; ind < m_fixedDimsVect.size(); ind++)
    {
      int idir = m_fixedDimsVect[ind];
      returnBox.setRange(idir, m_fixedPt[ind]);
    }
  return returnBox;
}


IntVectSet
MultiBlockUtil::extraBlockGhosts(const Box& a_baseBox,
                                 int a_ghostLayer,
                                 int a_baseBlockNum) const
{
  CH_assert(m_isDefined);
  const Vector<Box>& mappingBlocksAll = m_coordSysPtr->mappingBlocks();
  const Box& baseBlockBox = mappingBlocksAll[a_baseBlockNum];

  // Without m_fixedDims: Box grownBox = grow(a_baseBox, a_ghostLayer);
  // Grow box in the interpolating dimensions only.
  IntVect ghostLayerVect = IntVect::Zero;
  for (int ind = 0; ind < m_interpDimsVect.size(); ind++)
    {
      int idir = m_interpDimsVect[ind];
      ghostLayerVect[idir] = a_ghostLayer;
    }
  Box grownBox = grow(a_baseBox, ghostLayerVect);
  
  IntVectSet ivsReturn = IntVectSet(); // empty
  if ( !baseBlockBox.contains(grownBox) ) // otherwise, empty
    {
      // We'll set ivs to extra-block ghost cells of baseBox.
      // Start with ivs being the full grownBox.
      DenseIntVectSet ivs = DenseIntVectSet(grownBox, true);
      // Note baseBlockBox includes all of baseBox,
      // so after this, ivs contains ghost cells only.
      ivs -= baseBlockBox;
      if ( !ivs.isEmpty() )
        {
          // Now remove ghost cells that are outside domain.
          ivs.recalcMinBox();
          // WAS for (int idir = 0; idir < SpaceDim; idir++)
          for (int ind = 0; ind < m_interpDimsVect.size(); ind++)
            {
              int idir = m_interpDimsVect[ind];
              // Check idir Lo face.
              int endBlockLo = baseBlockBox.smallEnd(idir);
              if (grownBox.smallEnd(idir) < endBlockLo)
                { // Remove ghost cells from idir Lo face.
                  if (m_boundaries[a_baseBlockNum][idir].isDomainBoundary())
                    { // Remove ghost cells beyond this face.
                      Box grownBoxFace(grownBox);
                      grownBoxFace.setBig(idir, endBlockLo-1);
                      ivs -= grownBoxFace;
                    }
                }
              // Check idir Hi face.
              int endBlockHi = baseBlockBox.bigEnd(idir);
              if (grownBox.bigEnd(idir) > endBlockHi)
                { // Remove ghost cells from idir Hi face.
                  if (m_boundaries[a_baseBlockNum][idir + SpaceDim].isDomainBoundary())
                    { // Remove ghost cells beyond this face.
                      Box grownBoxFace(grownBox);
                      grownBoxFace.setSmall(idir, endBlockHi+1);
                      ivs -= grownBoxFace;
                    }
                }
            }
          ivs.recalcMinBox();
          // Now ivs is what we want it to be.
          // Find the valid block and valid cell of each ghost
          // cell in ivs.
          // I have to convert from DenseIntVectSet.
          ivsReturn = IntVectSet(ivs);
        }
    }
  return ivsReturn;
}

bool
MultiBlockUtil::allGridsHaveFixedPt(const BoxLayout& a_layout)
{
  if (m_fixedDimsVect.size() == 0)
    {
      return true;
    }
  else
    {
      for (LayoutIterator lit = a_layout.layoutIterator(); lit.ok(); ++lit)
        {
          const Box& bx = a_layout[lit];
          for (int ind = 0; ind < m_fixedDimsVect.size(); ind++)
            {
              int idir = m_fixedDimsVect[ind];
              int val = m_fixedPt[ind];
              if ( ! ( ( bx.smallEnd(idir) <= val ) &&
                       ( val <= bx.bigEnd(idir) ) ) )
                { // val is not in range of bx[idir]
                  return false;
                }
            }
        }
      // every box passed the test in every fixed dimension
      return true;
    }
}


void
MultiBlockUtil::getCollapsedLayout(BoxLayout& a_layoutCollapsed,
                                   const BoxLayout& a_layoutFull)
{
  if (m_fixedDimsVect.size() == 0)
    {
      a_layoutCollapsed = a_layoutFull;
    }
  else
    {
      BoxCollapser collapser(m_fixedDims);
      a_layoutCollapsed.deepCopy(a_layoutFull);
      a_layoutCollapsed.transform(collapser);
      a_layoutCollapsed.closeNoSort();
    }
}

void
MultiBlockUtil::getCollapsedLayout(DisjointBoxLayout& a_layoutCollapsed,
                                   const DisjointBoxLayout& a_layoutFull)
{
  if (m_fixedDimsVect.size() == 0)
    {
      a_layoutCollapsed = a_layoutFull;
    }
  else
    {
      BoxCollapser collapser(m_fixedDims);
      a_layoutCollapsed.deepCopy(a_layoutFull);
      a_layoutCollapsed.transform(collapser);
      a_layoutCollapsed.closeNoSort();
    }
}

void
MultiBlockUtil::getFixedOffLayout(BoxLayout& a_layoutFixedOff,
                                  const BoxLayout& a_layoutFull)
{
  BoxFixedOff fixer(m_fixedDims);
  a_layoutFixedOff.deepCopy(a_layoutFull);
  a_layoutFixedOff.transform(fixer);
  a_layoutFixedOff.closeNoSort();
}

void
MultiBlockUtil::order2grad(LevelData<FArrayBox>&        a_gradData,
                           const LevelData<FArrayBox>&  a_data)
{
  CH_assert(m_isDefined);
  // This function fills in a_gradData on valid cells ONLY.
  const DisjointBoxLayout& layout = a_data.disjointBoxLayout();
  int ncomp = a_data.nComp();
  LevelData<FArrayBox> dataGhosted(layout, ncomp, m_interpUnit);
  DataIterator dit = a_data.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& baseBox = layout[dit];

      const FArrayBox& dataFab = a_data[dit];
      FArrayBox& dataGhostedFab = dataGhosted[dit];

      dataGhostedFab.copy(dataFab);
      // In each interpolating dimension, fill a layer of ghost cells
      // with extrapolant.
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if ( !m_fixedDims.contains(idir) )
            {
              secondOrderCellExtrap(dataGhostedFab, baseBox, idir);
            }
        }
    }
  // Where it exists, valid data from neighboring boxes overwrites
  // the extrapolated data in dataGhosted.
  // If we have data covering the whole domain, then all extrapolated
  // data will be overwritten.  But the usefulness here comes when
  // we do not have data covering the whole domain.  It works out
  // that we apply the second-order one-sided finite-difference
  // derivative formula on the edge of where we have data.
  dataGhosted.exchange();

  // Now find gradient from dataGhosted.
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& baseBox = layout[dit];
      ProblemDomain blockDomain = m_coordSysPtr->blockDomainOfBox(baseBox);

      const FArrayBox& dataGhostedFab = dataGhosted[dit];
      FArrayBox& gradFab = a_gradData[dit];

      // Don't use any of dataGhostedFab from outside blockDomain.
      Box dataBox(dataGhostedFab.box());
      dataBox &= blockDomain;
      
      // This function is in MBMiscUtil.
      order2gradient(gradFab, baseBox,
                     dataGhostedFab, dataBox,
                     m_interpDimsVect);
    }

  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused++;
}

#include "NamespaceFooter.H"
