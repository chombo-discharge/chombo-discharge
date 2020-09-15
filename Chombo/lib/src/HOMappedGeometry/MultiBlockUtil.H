#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MULTIBLOCKUTIL_H_
#define _MULTIBLOCKUTIL_H_

#include "MultiBlockCoordSys.H"
#include "IndicesTransformation.H"
#include "IVSFAB.H"
#include "MBStencil.H"
#include "MBVectorStencil.H"

#include "NamespaceHeader.H"

/// The class MultiBlockUtil contains some functions used by MultiBlockLevelGeom, MultiBlockLevelExchange, and MultiBlockLevelCopier. It does NOT store any layouts.
class MultiBlockUtil
{
public:

  /// default constructor
  MultiBlockUtil();

  /// destructor
  virtual ~MultiBlockUtil();

  /// full constructor
  MultiBlockUtil(const MultiBlockCoordSys*  a_coordSysPtr,
                 Interval                   a_fixedDims = Interval(),
                 Vector<int>                a_fixedPt = Vector<int>() );

  virtual void undefine();

  virtual void define(const MultiBlockCoordSys*  a_coordSysPtr,
                      Interval                   a_fixedDims = Interval(),
                      Vector<int>                a_fixedPt = Vector<int>() );

  /// the offset in [-1:1]^D of cell a_cell from the box of block a_block.
  IntVect whichOffset(int             a_block,
                      const IntVect&  a_cell) const;

  /// the transformation from block a_dstBlock to a_srcBlock, taking offset a_offset from a_dstBlock
  const IndicesTransformation& blockTransform(int a_srcBlock,
                                              int a_dstBlock,
                                              const IntVect& a_offset) const;

  /// transformations to block a_srcBlock, taking offset a_offset from all destination blocks
  const Vector<IndicesTransformation>& transformsSrcOff(int a_srcBlock,
                                                        const IntVect& a_offset) const;

  /// add IntVects to an IntVectSet in all directions, even across block boundaries.
  void growIVS(IntVectSet& a_ivs,
               int a_bufferSize) const;

  /// return list of all cells that share a vertex with a given cell
  /** Given a_baseCell in block a_baseBlock, return a Vector of cells
      with indices a_neighborCells in corresponding blocks a_neighborBlocks.
  */
  void commonVertexCells(Vector<IntVect>& a_neighborCells,
                         Vector<int>& a_neighborBlocks,
                         const IntVect& a_baseCell,
                         int a_baseBlock) const;

  void commonVertexStencilElements(Vector<MBStencilElement>& a_neighbors,
                                   const MBStencilElement& a_base) const;

  /// return shifted cell
  /**
     Starting with valid cell with indices a_origCell in block a_origBlock;
     shift in direction a_dir by a_shift;
     end up in valid cell with indices a_shiftedCell in block a_shiftedBlock.
     Return false if go beyond physical boundary.
   */
  bool validCellShift(IntVect&        a_shiftedCell,
                      int&            a_shiftedBlock,
                      const IntVect&  a_origCell,
                      int             a_origBlock,
                      int             a_dir,
                      int             a_shift) const;

  bool validCellShift(MBStencilElement& a_shiftedElt,
                      const MBStencilElement& a_origElt,
                      int             a_dir,
                      int             a_shift) const;

  /// return index, block, and center in mapped coordinates of valid cell containing each ghost cell
  /**
     For each IntVect thisGhostCell in a_ghostCellsIVS, which is
     an extra-block ghost cell of block a_baseBlockNum,
     find the valid cell that contains its center.
     Then set a_validIndices(thisGhostCell, 0) to its indices,
     set a_validBlock(thisGhostCell, 0) to its block number,
     and set a_validMappedCenter(thisGhostCell, 0) to the center of
     thisGhostCell in the mapped coordinates of its valid block.
   */
  void getValid(IVSFAB<IntVect>& a_validIndices,
                IVSFAB<int>& a_validBlock,
                IVSFAB<RealVect>& a_validMappedCenter,
                const IntVectSet& a_ghostCellsIVS,
                int a_baseBlockNum) const;

  void getValid(BaseFab<IntVect>& a_validIndices,
                BaseFab<int>& a_validBlock,
                BaseFab<RealVect>& a_validMappedCenter,
                const Box& a_bx,
                int a_baseBlockNum) const;

  /// find the stencil cells and source blocks for all cells in a_ghostCellsIVS, and also the minimum Box covering all such cells for each source block
  /**
     In a_stencil, return the valid cells and source block numbers of the
     stencil of radius a_radius for every ghost cell in a_ghostCellsIVS.
     In a_minBox[srcBlock], return the minimum Box that covers all
     the cells of a_stencils that are in block srcBlock.

     For each ghost cell thisGhostCell in a_ghostCellsIVS,
     the center of thisGhostCell lies in block a_validBlock(thisGhost, 0)
     within the valid cell with indices a_validIndices(thisGhostCell, 0)
     and within that block, the ghost cell center has mapped coordinates
     a_validMappedCenter(thisGhostCell, 0).
   */
  void getStencilCells(IVSFAB<MBStencil>& a_stencils,
                       Vector<Box>& a_minBox,
                       const IntVectSet& a_ghostCellsIVS,
                       const IVSFAB<IntVect>& a_validIndices,
                       const IVSFAB<int>& a_validBlock,
                       const IVSFAB<RealVect>& a_validMappedCenter,
                       int a_order,
                       int a_radius) const;

  void getVectorStencilCells(IVSFAB<MBVectorStencil>& a_vectorstencils,
                             const IVSFAB<MBStencil>& a_stencils,
                             const IntVectSet& a_ghostCellsIVS) const;

  void getVectorTransformations(Vector< BaseFab<VectorTransformation>* >& a_vectorTransformations,
                                const IVSFAB<MBVectorStencil>& a_vectorstencilsFab,
                                const IndexType& a_type,
                                const IntVectSet& a_ghostCellsIVS,
                                int a_ghostBlockNum) const;

  void getWeights(IVSFAB<MBStencil>& a_stencilsFab,
                  const IndexType& a_type,
                  const IntVectSet& a_ghostCellsIVS,
                  const IVSFAB<IntVect>& a_validIndices,
                  const IVSFAB<int>& a_validBlock,
                  const IVSFAB<RealVect>& a_validMappedCenter,
                  const Vector<IntVect>& a_exponents,
                  int a_ghostBlockNum) const;

  void getVectorWeights(IVSFAB<MBVectorStencil>& a_vectorstencilsFab,
                        const Vector< BaseFab<VectorTransformation>* >& a_vectorStencilTransformations,
                        const IndexType& a_type,
                        const IntVectSet& a_ghostCellsIVS,
                        const IVSFAB<IntVect>& a_validIndices,
                        const IVSFAB<int>& a_validBlock,
                        const IVSFAB<RealVect>& a_validMappedCenter,
                        const Vector<IntVect>& a_exponents,
                        int a_ghostBlockNum) const;

  /// copy from a_blockStencilsFab to a_stencilsFab on a_ghostCellsIVS, where the cells in the stencil of a_blockStencilsFab are m_fixedPt in m_fixedDimsVect.
  void copyStencilsFromBlock(IVSFAB<MBStencil>& a_stencilsFab,
                             const IntVectSet& a_ghostCellsIVS,
                             const IVSFAB<MBStencil>& a_blockStencilsFab) const;

  void copyVectorStencilsFromBlock(IVSFAB<MBVectorStencil>& a_vectorstencilsFab,
                                   const IntVectSet& a_ghostCellsIVS,
                                   const IVSFAB<MBVectorStencil>& a_blockVectorStencilsFab) const;

  Box boxFixed(const Box& a_bx);

  IntVectSet extraBlockGhosts(
                              const Box& a_baseBox,
                              int a_ghostLayer,
                              int a_baseBlockNum) const;

  bool allGridsHaveFixedPt(const BoxLayout& a_layout);

  void getCollapsedLayout(BoxLayout& a_layoutCollapsed,
                          const BoxLayout& a_layoutFull);

  void getCollapsedLayout(DisjointBoxLayout& a_layoutCollapsed,
                          const DisjointBoxLayout& a_layoutFull);

  void getFixedOffLayout(BoxLayout& a_layoutFixedOff,
                         const BoxLayout& a_layoutFull);

  /// returns a_gradData = gradient(a_data) on valid cells, taking 1-sided differences next to block boundaries
  void order2grad(LevelData<FArrayBox>&        a_gradData,
                  const LevelData<FArrayBox>&  a_data);

protected:

  /// return list of all aliases of a vertex
  /** Given a_baseVertex, a vertex of block a_baseBlock,
      return a Vector of vertices with indices a_aliasVertices 
      in corresponding blocks a_aliasBlocks.
  */
  void vertexAliases(Vector<IntVect>& a_aliasVertices,
                     Vector<int>& a_aliasBlocks,
                     const IntVect& a_baseVertex,
                     int a_baseBlock) const;

  void vertexStencilElementAliases(Vector<MBStencilElement>& a_aliases,
                                   const MBStencilElement& a_base) const;

  /// return the set of all valid cells adjacent to a given set of vertices
  /** Given list of vertices a_vertexIndices in blocks a_vertexBlocks,
      return a Vector of cell indices a_cellIndices
      in corresponding blocks a_cellBlocks.
  */
  void validCellsFromVertices(Vector<IntVect>& a_cellIndices,
                              Vector<int>& a_cellBlocks,
                              const Vector<IntVect>& a_vertexIndices,
                              const Vector<int>& a_vertexBlocks) const;

  void validStencilElementsFromVertices(Vector<MBStencilElement>& a_stencilElements,
                                        const Vector<MBStencilElement>& a_vertices) const;

  /// return <((Xvec() - Xvec(a_commonMappedBasePoint))/a_avgDistance)^pvec>_(a_cell), for all pvec
  /**
     Sets a_powers[pvec] = <((Xvec() - Xvec(a_commonMappedBasePoint)/a_avgDistance))^pvec>_(a_cell)
     for each pvec such that sum(|pvec|) <= m_degree
     where
     a_cell is in index space of a_thisBlockNum (may or may not be valid here),
     XVec() is in mapped space of a_commonBlockNum and varies over a_cell,
     a_commonMappedBasePoint is a point in mapped space of a_commonBlockNum .

     Example in 2D with order 4, in this sequence:
     0: 0 0
     1: 1 0
     2: 2 0
     3: 3 0
     4: 0 1
     5: 1 1
     6: 2 1
     7: 0 2
     8: 1 2
     9: 0 3
   */
  void displacementPowers(Vector<Real>&      a_powers,
                          const IntVect&     a_cell,
                          const RealVect&    a_commonMappedBasePoint,
                          Real               a_avgDistance,
                          int                a_thisBlockNum,
                          int                a_commonBlockNum,
                          const Vector<IntVect>&   a_exponents) const;

  void displacementPowersTransformed(Vector<Real>&      a_powers,
                                     const BaseFab<VectorTransformation>& a_vectorStencilTransformations,
                                     const IntVect&     a_cell,
                                     const RealVect&    a_commonMappedBasePoint,
                                     Real               a_avgDistance,
                                     int                a_thisBlockNum,
                                     int                a_commonBlockNum,
                                     int                a_ghostBlockNum,
                                     const Vector<IntVect>&   a_exponents) const;

  void getStencilCells(MBStencil& a_stencil,
                       Vector<Box>& a_minBox,
                       const IntVect& a_iv,
                       const IntVect& a_validIndices,
                       int a_validBlockNum,
                       const RealVect& a_validMappedCenter,
                       int a_order,
                       int a_radius) const;

  void getWeights(MBStencil& a_stencil,
                  const IndexType& a_type,
                  const IntVect& a_iv,
                  const IntVect& a_validIndices,
                  int a_validBlockNum,
                  const RealVect& a_validMappedCenter,
                  const Vector<IntVect>& a_exponents,
                  int a_ghostBlockNum) const;

  void getVectorWeights(MBVectorStencil& a_vectorstencil,
                        const Vector< BaseFab<VectorTransformation>* >& a_vectorStencilTransformations,
                        const IndexType& a_type,
                        const IntVect& a_iv,
                        const IntVect& a_validIndices,
                        int a_validBlockNum,
                        const RealVect& a_validMappedCenter,
                        const Vector<IntVect>& a_exponents,
                        int a_ghostBlockNum) const;

  void copyStencilsFromBlock(MBStencil& a_stencil,
                             const IntVect& a_iv,
                             const MBStencil& a_blockStencil) const;

  void copyVectorStencilsFromBlock(MBVectorStencil& a_vectorstencil,
                                   const IntVect& a_iv,
                                   const MBVectorStencil& a_blockVectorStencil) const;

  /// is defined?
  bool m_isDefined;

  /// contains coordinate system
  MultiBlockCoordSys* m_coordSysPtr;

  /// all mapping blocks, from m_coordSysPtr
  Vector<Box> m_mappingBlocks;

  /// Number of mapping blocks
  int m_nblocks;

  /// all block boundaries; from m_coordSysPtr
  Vector< Tuple<BlockBoundary, 2*SpaceDim> > m_boundaries;

  Box m_offsetAllBox;

  /// transformations from destination blocks to source blocks. Has m_nblocks components and lives on m_offsetAllBox.
  /**
     IndicesTransformation (*m_transformsAll(ivOff, srcBlock)[dstBlock]
     is the transformation from dstBlock to srcBlock
     in the chunk of the halo offset by IntVect ivOff from dstBlock.
  */
  // Vector< BaseFab<IndicesTransformation>* > m_transformsAll;
  BaseFab< Vector<IndicesTransformation>* > m_transformsAll;

  /// dimensions that are fixed.  Default is empty.
  Interval m_fixedDims;

  /// dimensions that are fixed.  Default is empty.
  Vector<int> m_fixedDimsVect;

  /// coordinates in dimensions of m_fixedDims.  Default is empty.
  Vector<int> m_fixedPt;

  /// dimensions that will be interpolated, which is all dimensions that are not in m_fixedDims.  Default is 0:SpaceDim-1.
  Vector<int> m_interpDimsVect;

  /// 1 in interpolated dimensions, 0 in fixed dimensions
  IntVect m_interpUnit;

  /// these are used in the the (2*SpaceDim+1)-point 4th-order convolution stencil
  int m_convSize;
  Vector<IntVect> m_convBaseCells;
  Vector<Real> m_convWeight;
};

#include "NamespaceFooter.H"

#endif // include guard
