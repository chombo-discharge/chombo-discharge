#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MULTIBLOCKLEVELGEOM_H_
#define _MULTIBLOCKLEVELGEOM_H_

#include "MultiBlockCoordSys.H"
#include "MultiBlockUtil.H"
#include "IVSFAB.H"
#include "LevelData.H"

#include "NamespaceHeader.H"

/// The class MultiBlockLevelGeom stores topological and geometric information
/**
   MultiBlockLevelGeom stores topological and geometric information
   about a particular layout with a particular MultiBlockCoordSys.

   Implementations of specific multiblock coordinate systems will be
   provided by derived classes.
*/
class MultiBlockLevelGeom
{

public:

  /// default constructor
  MultiBlockLevelGeom()
  {
    m_isDefined = false;
  }

  /// destructor
  ~MultiBlockLevelGeom();

  /// full constructor
  MultiBlockLevelGeom(const MultiBlockCoordSys*  a_coordSysPtr,
                      const DisjointBoxLayout&   a_grids,
                      int                        a_ghosts,
                      Interval                   a_fixedDims = Interval(),
                      Vector<int>                a_fixedPt = Vector<int>());

  void undefine();

  /// define knowing only ghosts required for the solution update procedure. DEPRECATED.
  void define(const MultiBlockCoordSys* a_coordSysPtr,
              const DisjointBoxLayout&  a_grids,
              const int                 a_solUpdateGhosts,
              const int                 a_spaceOrder,
              Interval                  a_fixedDims = Interval(),
              Vector<int>               a_fixedPt = Vector<int>());


  /// define with explicit specification of number of multiblock ghosts
  void define(const MultiBlockCoordSys*  a_coordSysPtr,
              const DisjointBoxLayout&   a_grids,
              int                        a_ghosts,
              Interval                   a_fixedDims = Interval(),
              Vector<int>                a_fixedPt = Vector<int>());

  /// accessor
  int ghosts() const
  { return m_ghosts; }

  /// accessor
  const DisjointBoxLayout& grids() const
  { return m_grids; }

  /// accessor
  const DisjointBoxLayout& gridsFull() const
  { return m_gridsFull; }

  /// accessor
  MultiBlockCoordSys* coordSysPtr() const
  { return m_coordSysPtr; }

  /// accessor
  const LayoutData<int>& block() const
  { return m_block; }

  /// accessor
  const LayoutData< IVSFAB<IntVect>* >& validIndices() const
  { return m_validIndices; }

  /// accessor
  const LayoutData< IVSFAB<int>* >& validBlock() const
  { return m_validBlock; }

  /// accessor
  const LayoutData< IVSFAB<RealVect>* >& validMappedCenter() const
  { return m_validMappedCenter; }

  const Interval& fixedDims() const
  { return m_fixedDims; }

  const Vector<int>& fixedDimsVect() const
  { return m_fixedDimsVect; }

  const Vector<int>& interpDimsVect() const
  { return m_interpDimsVect; }

  const Vector<int>& fixedPt() const
  { return m_fixedPt; }

  const BoxLayout& gridsFixedOff() const
  { return m_gridsFixedOff; }

  bool allGridsHaveFixedPt() const
  { return m_allGridsHaveFixedPt; }

  MultiBlockUtil* mbUtil() const
  { return m_mbUtil; }

  /// return indices of destination ghost cells for a particular width
  /**
    In a_ghostCells, return ghost cells of m_grids out to a_ghostLayer,
    which is width of destination ghost cell layer.
    This function requires
    m_grids, m_coordSysPtr, m_block, and m_boundaries, but
    NOT anything else like m_validIndices, m_validBlock, m_validMappedCenter.
  */
  void extraBlockGhosts(LayoutData< IntVectSet >& a_ghostCells,
                        int a_ghostLayer) const;

  bool isDefined() const
  {
    return m_isDefined;
  }
  
protected:

  /// is this object defined?
  /// covers m_coordSysPtr, m_grids, m_ghosts, m_boundaries, m_block, m_ghostCells, m_validIndices, m_validBlock, m_validMappedCenter
  bool m_isDefined;

  /// coordinate system (given in constructor)
  MultiBlockCoordSys* m_coordSysPtr;

  /// original underlying grids
  DisjointBoxLayout m_gridsFull;

  /// underlying grids (given in constructor) varying in m_interpDimsVect but set to m_fixedPt in m_fixedDimsVect
  DisjointBoxLayout m_grids;

  /// displacement of original underlying grids from m_grids: 0 in m_interpDimsVect but varying in m_fixedDimsVect
  BoxLayout m_gridsFixedOff;

  /// whether all original underlying grids contain m_fixedPt in m_fixedDimsVect
  bool m_allGridsHaveFixedPt;

  /// width of ghost layer around m_grids, from which data will be taken (given in constructor)
  int m_ghosts;

  /// dimensions that are fixed.  Default is empty.
  Interval m_fixedDims;

  /// dimensions that are fixed.  Default is empty.
  Vector<int> m_fixedDimsVect;

  /// coordinates in dimensions of m_fixedDims.  Default is empty.
  Vector<int> m_fixedPt;

  /// dimensions that will be interpolated, which is all dimensions that are not in m_fixedDims.  Default is 0:SpaceDim-1.
  Vector<int> m_interpDimsVect;

  /// all block boundaries; from m_coordSysPtr
  Vector< Tuple<BlockBoundary, 2*SpaceDim> > m_boundaries;

  /// block containing each box of m_grids
  LayoutData<int> m_block;

  /// extra-block ghost cells from which data will be taken.  Lives on m_grids.
  LayoutData< IntVectSet > m_ghostCells;

  /// coordinate indices within the valid block of the cell in which the center of each extra-block ghost cell falls.  Lives on m_ghostCells, has 1 component.
  LayoutData< IVSFAB<IntVect>* > m_validIndices;

  /// index of the valid block in which the center of each extra-block ghost cell falls.  Lives on m_ghostCells, has 1 component.
  LayoutData< IVSFAB<int>* > m_validBlock;

  /// mapped coordinates, in valid block's index space, of center of each extra-block ghost cell.  Lives on m_ghostCells, has 1 component.
  LayoutData< IVSFAB<RealVect>* > m_validMappedCenter;

  MultiBlockUtil* m_mbUtil;
};

#include "NamespaceFooter.H"
#endif
