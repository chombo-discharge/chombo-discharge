#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MBAggLevelExchange__H_
#define _MBAggLevelExchange__H_

#include "MultiBlockLevelGeom.H"
#include "IVSFAB.H"
#include "VectorTransformation.H"
#include "MBStencil.H"
#include "MBAggStencil.H"
#include "NamespaceHeader.H"

/// The class MultiBlockLevelExchange stores topological and geometric information
/**
   MultiBlockLevelExchange stores topological and geometric information
   about a particular layout with a particular MultiBlockCoordSys.

   Implementations of specific multiblock coordinate systems will be
   provided by derived classes.
*/
class MBAggLevelExchange
{

public:

  /// default constructor
  MBAggLevelExchange()
  { m_isDefined = false;
    m_isVectorDefined = false; }

  /// destructor
  virtual ~MBAggLevelExchange();

  /// full constructor
  MBAggLevelExchange(const MultiBlockLevelGeom*  a_geomPtr,
                     int                         a_ghosts,
                     int                         a_order);

  virtual void undefine();

  virtual void define(const MultiBlockLevelGeom*  a_geomPtr,
                      int                         a_ghosts,
                      int                         a_order);

  /// must be in the derived class, else error when calling.
  virtual void defineVector();

  /// fill extra-block ghost cells of a_data at components in a_intvl
  virtual void interpGhosts(LevelData<FArrayBox>&  a_data,
                            const Interval&        a_intvl) const;

  /// fill extra-block ghost cells of a_data at all components
  virtual void interpGhosts(LevelData<FArrayBox>&  a_data) const;

  /// fill extra-block ghost cells of vector a_data at all components
  virtual void interpGhostsVector(LevelData<FArrayBox>&  a_data) const;

  /// fill extra-block ghost cells of vector a_data at components in a_intvl, of length SpaceDim.  Must be in the derived class, else error when calling.
  virtual void interpGhostsVector(LevelData<FArrayBox>&  a_data,
                                  const Interval&        a_intvl) const;

  /// fill extra-block ghost cells of all components of a_data at components, where a_vecIntvl is interval of components of a vector, and the other components are scalars
  virtual void interpGhostsAllWithVector(LevelData<FArrayBox>&  a_data,
                                         const Interval&        a_vecIntvl) const;

  /// accessor
  const LayoutData< IntVectSet >& ghostCells() const
  { CH_assert(isDefined());
    return m_ghostCells; }

  /// accessor

  /// accessor
  /**
     IVSFAB<MBStencil> *m_stencils[dit] lives on m_ghostCells[dit]
     and has 1 component.
     For each IntVect thisGhostCell in IntVectSet m_ghostCells[dit]
     (i.e., thisGhostCell an extra-block ghost cell of m_grids[dit]),
     MBStencil *m_stencils[dit](thisGhostCell, 0) contains the
     interpolation stencil for thisGhostCell.
  */
  const LayoutData< RefCountedPtr< MBAggStencil > >& stencils() const
  { CH_assert(isDefined());
    return m_stencils; }

  bool isDefined() const
  {
    return m_isDefined;
  }

protected:

  /// remove from m_ghostCells all cells that do not have complete stencils.
  virtual void removeNoValidSource(const LayoutData< RefCountedPtr< IVSFAB<MBStencil> > >& a_slowStencils);


  /// is defined?
  bool m_isDefined;

  /// is defined for vector transformations?
  bool m_isVectorDefined;

  /// either IndexType::TheCellType() for averaged, or IndexType::TheNodeType() for centered
  IndexType m_type;

  /// underlying grids
  DisjointBoxLayout m_grids;

  /// width of ghost layer around m_grids, at which data will be filled in
  int m_ghosts;

  /// order of approximation wanted
  int m_order;

  /// radius of neighborhood of data used for interpolation to ghost cell
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

  /// original underlying grids
  DisjointBoxLayout m_gridsFull;

  /// Number of mapping blocks
  int m_nblocks;

  /// extra-block ghost cells of m_grids that need to be filled in
  LayoutData< IntVectSet > m_ghostCells;

  /// interval of dimensions that are fixed.  Default empty.
  Interval m_fixedDims;

  /// dimensions that are fixed.  Default empty.
  Vector<int> m_fixedDimsVect;

  /// dimensions that are not fixed.  Default 0:SpaceDim-1.
  Vector<int> m_interpDimsVect;

  /// fixed values of dimensions that are fixed.  Default empty.
  Vector<int> m_fixedPt;

  /// displacement of original underlying grids from m_grids: 0 in m_interpDimsVect but varying in m_fixedDimsVect
  BoxLayout m_gridsFixedOff;

  /// whether all original underlying grids contain m_fixedPt in m_fixedDimsVect
  bool m_allGridsHaveFixedPt;


  /// interpolation stencils
  /**
     For each dit, IVSFAB<MBStencil> *m_stencils[dit]
     lives on IntVectSet m_ghostCells[dit] and has 1 component,
     containing the stencil for each extra-block ghost cell of m_grids[dit].
  */
  LayoutData< RefCountedPtr< MBAggStencil > > m_stencils;

  /// for each source block, min Box containing valid cells
  /**
     For each srcBlock and dit,
     Box (*(m_validLayout[srcBlock]))[dit]
     is the minimum Box containing all valid cells that are needed
     to fill in the ghost cells of m_ghostCells[dit].
  */
  Vector<BoxLayout*> m_validLayout;

  Vector<BoxLayout*> m_validFullLayout;

  /// Copier *m_copiers[srcBlock] for copying from m_grids to BoxLayout *m_validLayout[srcBlock]
  Vector<Copier*> m_copiers;

  Vector<Copier*> m_copiersFull;

  MultiBlockUtil* m_mbUtil;
};

#include "NamespaceFooter.H"
#endif
