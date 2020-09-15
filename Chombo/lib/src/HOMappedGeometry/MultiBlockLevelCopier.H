#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MULTIBLOCKLEVELCOPIER_H_
#define _MULTIBLOCKLEVELCOPIER_H_

#include "MultiBlockLevelGeom.H"
#include "IVSFAB.H"
#include "MultiBlockUtil.H"
#include "MBStencil.H"

#include "NamespaceHeader.H"

/// The class MultiBlockLevelCopier stores topological and geometric information
/**
   MultiBlockLevelCopier stores topological and geometric information
   about a particular layout with a particular MultiBlockCoordSys,
   for copying to another layout.

   Implementations of specific multiblock coordinate systems will be
   provided by derived classes.
*/
class MultiBlockLevelCopier
{

public:

  /// default constructor
  MultiBlockLevelCopier()
  { m_isDefined = false;
    m_isVectorDefined = false; }

  /// destructor
  virtual ~MultiBlockLevelCopier();

  /// full constructor
  MultiBlockLevelCopier(const MultiBlockLevelGeom*  a_geomPtr,
                        const BoxLayout&            a_dstLayout,
                        int                         a_order);

  MultiBlockLevelCopier(const MultiBlockLevelGeom*  a_geomPtr,
                        const DisjointBoxLayout&    a_dstDisjointLayout,
                        const IntVect&              a_ghostVect,
                        int                         a_order);

  virtual void undefine();

  virtual void define(const MultiBlockLevelGeom*  a_geomPtr,
                      const BoxLayout&            a_dstLayout,
                      int                         a_order);

  virtual void define(const MultiBlockLevelGeom*  a_geomPtr,
                      const DisjointBoxLayout&    a_dstDisjointLayout,
                      const IntVect&              a_ghostVect,
                      int                         a_order);

  /// must be in the derived class, else error when calling.
  virtual void defineVector();

  /// copy from a_src to extra-block cells of a_dst at components in a_intvl
  virtual void copyTo(const LevelData<FArrayBox>&  a_src,
                      BoxLayoutData<FArrayBox>&    a_dst,
                      const Interval&        a_intvl) const;

  /// copy from a_src to extra-block cells of a_dst at all components
  virtual void copyTo(const LevelData<FArrayBox>&  a_src,
                      BoxLayoutData<FArrayBox>&    a_dst) const;

  /// copy from a_src to extra-block cells of vector a_dst at all components
  virtual void copyToVector(const LevelData<FArrayBox>&  a_src,
                            BoxLayoutData<FArrayBox>&    a_dst) const;

  /// copy from a_src to extra-block cells of vector a_dst at components in a_intvl, of length SpaceDim.  Must be in the derived class, else error when calling.
  virtual void copyToVector(const LevelData<FArrayBox>&  a_src,
                            BoxLayoutData<FArrayBox>&    a_dst,
                            const Interval&        a_intvl) const;

  /// copy from a_src to extra-block cells of a_dst at all components where a_vecIntvl is interval of components of a vector, and the other components are scalars
  virtual void copyToAllWithVector(const LevelData<FArrayBox>&  a_src,
                                   BoxLayoutData<FArrayBox>&    a_dst,
                                   const Interval&         a_vecIntvl) const;


  /// accessor
  const LayoutData< IntVectSet >& extraCells() const
  { CH_assert(isDefined());
    return m_extraCellsIVS; }

  /// accessor
  /**
     IVSFAB<MBStencil> *m_stencils[dit] lives on m_ghostCells[dit]
     and has 1 component.
     For each IntVect thisGhostCell in IntVectSet m_ghostCells[dit]
     (i.e., thisGhostCell an extra-block ghost cell of m_grids[dit]),
     MBStencil *m_stencils[dit](thisGhostCell, 0) contains the
     interpolation stencil for thisGhostCell.
   */
  const LayoutData< RefCountedPtr< IVSFAB<MBStencil> > >& stencils() const
  { CH_assert(isDefined());
    return m_stencils; }

  bool isDefined() const
  {
    return m_isDefined;
  }

protected:

  /// is defined?
  bool m_isDefined;

  /// is defined for vector transformations?
  bool m_isVectorDefined;

  /// either IndexType::TheCellType() for averaged, or IndexType::TheNodeType() for centered
  IndexType m_type;

  /// underlying source grids
  DisjointBoxLayout m_srcLayout;

  DisjointBoxLayout m_srcFullLayout;

  /// underlying destination grids
  BoxLayout m_dstLayout;

  /// displacement of original destination layout from m_dstLayout: 0 in m_interpDimsVect but varying in m_fixedDimsVect
  BoxLayout m_dstFixedOff;

  /// whether all original underlying grids contain m_fixedPt in m_fixedDimsVect
  bool m_allGridsHaveFixedPt;

  /// block number of base box of each box in m_dstLayout
  LayoutData<int> m_blockNumbersDst;

  /// order of approximation wanted
  int m_order;

  /// radius of neighborhood of data used for interpolation to extra-block cell
  int m_radius;

  /// degree of Taylor polynomials
  int m_degree;

  /// [0:m_degree]^D
  Box m_degreeBox;

  /// powers in Taylor polynomial
  Vector<IntVect> m_powers;

  /// number of coefficients in Taylor polynomial
  int m_numTaylorCoeffs;

  /// contains topological and geometric information
  MultiBlockLevelGeom* m_geomPtr;

  /// contains coordinate systems; from m_geomPtr
  MultiBlockCoordSys* m_coordSysPtr;

  /// Number of mapping blocks
  int m_nblocks;

  /// extra-block cells of m_dstLayout that need to be filled in, in the form of IntVectSets
  /**
     For each dit, IntVectSet m_extraCellsIVS[dit]
     contains the cells of m_grids[dit] that\
     are not in block m_blockNumbersDst[dit].
   */
  LayoutData< IntVectSet > m_extraCellsIVS;

  /// index of the valid block in which the center of each extra-block ghost cell falls.  Lives on m_extraCellsIVS, has 1 component.
  LayoutData< IVSFAB<int>* > m_validBlock;

  /// index of the valid cell in which each extra-block ghost cell falls.  Lives on m_extraCellsIVS, has 1 component.
  LayoutData< IVSFAB<IntVect>* > m_validIndices;

  /// mapped coordinates, in valid block's index space, of center of each extra-block ghost cell.  Lives on m_extraCellsIVS, has 1 component.
  LayoutData< IVSFAB<RealVect>* > m_validMappedCenter;

  /// interval of dimensions that are fixed.  Default empty.
  Interval m_fixedDims;

  /// dimensions that are fixed.  Default empty.
  Vector<int> m_fixedDimsVect;

  /// dimensions that are not fixed.  Default 0:SpaceDim-1.
  Vector<int> m_interpDimsVect;

  /// fixed values of dimensions that are fixed.  Default empty.
  Vector<int> m_fixedPt;

  /// displacement from m_fixedPt of original underlying grids, 0 in m_interpDimsVect but varying in m_fixedDimsVect
  BoxLayout m_gridsFixedOff;

  /// interpolation stencils
  /**
     For each dit, IVSFAB<MBStencil> *m_stencils[dit]
     lives on IntVectSet m_extraCellsIVS[dit] and has 1 component,
     containing the stencil for each extra-block ghost cell of m_grids[dit].
  */
  LayoutData< RefCountedPtr< IVSFAB<MBStencil> > > m_stencils;

  /// vector interpolation stencils
  /**
     For each dit, IVSFAB<MBVectorStencil> *m_vectorstencils[dit]
     lives on IntVectSet m_extraCellsIVS[dit] and has 1 component,
     containing the vector stencil for each extra-block ghost cell of m_grids[dit].
  */
  LayoutData< RefCountedPtr< IVSFAB<MBVectorStencil> > > m_vectorstencils;

  /// for each source block, min Box containing valid cells
  /**
     For each srcBlock and dit,
     Box (*(m_validLayout[srcBlock]))[dit]
     is the minimum Box containing all valid cells that are needed
     to fill in the ghost cells of m_ghostCells[dit].
   */
  Vector<BoxLayout*> m_validLayout;

  Vector<BoxLayout*> m_validFullLayout;

  Vector< RefCountedPtr< LayoutData<Box> > > m_stencilCellsMinBox;

  Vector< RefCountedPtr< LayoutData<Box> > > m_stencilCellsFullMinBox;

  /// Copier *m_copiers[srcBlock] for copying from m_grids to BoxLayout *m_validLayout[srcBlock]
  Vector<Copier*> m_copiers;

  Vector<Copier*> m_copiersFull;

  MultiBlockUtil* m_mbUtil;
};

#include "NamespaceFooter.H"
#endif
