/*!
  @file jump_bc.H
  @brief Declaration of a multifluid boundary condition class for the poisson equation
  @author Robert Marskar
  @date Dec. 2017
*/

#ifndef _JUMP_BC_
#define _JUMP_BC_

#include "mfis.H"
#include "MFLevelGrid.H"
#include "MFInterfaceFAB.H"
#include "MFBaseIVFAB.H"

#include <Stencils.H>
#include <EBISBox.H>

/*!
  @brief Class for deriving the surface potential for multifluid poisson equations
*/
class jump_bc {
public:

  /*!
    @brief Empty constructor, must call define
  */
  jump_bc();

  /*!
    @brief Full constructor
  */
  jump_bc(const MFLevelGrid& a_mflg, const Real& a_dx, const int a_order, const LayoutData<IntVectSet>* a_cfivs);

  /*!
    @brief Destructor
  */
  virtual ~jump_bc();

  /*!
    @brief Switch between quadrant-based and non-quadrant based least squares stencils
  */
  static bool s_quadrant_based;

  /*!
    @brief Radius for least squares stencils
  */
  static int s_lsq_radius;

  /*!
    @brief Define function
  */
  virtual void define(const MFLevelGrid& a_mflg, const Real& a_dx, const int a_order, const LayoutData<IntVectSet>* a_cfivs);

  /*!
    @brief
  */
  virtual void match_bc(MFInterfaceFAB<Real>& m_phibc, const MFCellFAB& a_phi);


protected:

  RefCountedPtr<mfis> m_mfis;

  /*!
    @brief Level grid
  */
  MFLevelGrid m_mflg;

  /*!
    @brief Grid
  */
  DisjointBoxLayout m_grids;

  /*!
    @brief Problem domain
  */
  ProblemDomain m_domain;

  /*!
    @brief Stencil order
  */
  int m_order;

  /*!
    @brief Resolution at level
  */
  Real m_dx;

/*!
  @brief Is defined or not
*/
bool m_defined;

  /*!
    @brief Coarse-fine interface
  */
  const LayoutData<IntVectSet>* m_cfivs;

  /*!
    @brief Surface potential. Adding it on both sides for simplicity 
  */
  LayoutData<MFInterfaceFAB<Real> > m_soln;

  /*!
    @brief b-coefficient on each side
  */
  LayoutData<MFInterfaceFAB<Real> > m_bco;

  /*!
    @brief Weights on each side
  */
  LayoutData<MFInterfaceFAB<Real> > m_weights;

  /*!
    @brief Stencils on each side
  */
  LayoutData<MFInterfaceFAB<VoFStencil> > m_stencils;

  /*!
    @brief Build stencils for matching
  */
  virtual void build_stencils();

  /*!
    @brief Get a second order stencil 
  */
  virtual bool get_second_order_sten(Real&             a_weight,
				     VoFStencil&       a_stencil,
				     const VolIndex&   a_vof,
				     const EBISBox&    a_ebisbox,
				     const IntVectSet& a_cfivs);

  /*!
    @brief Get a first order stencil 
  */
  virtual void get_first_order_sten(Real&             a_weight,
				    VoFStencil&       a_stencil,
				    const VolIndex&   a_vof,
				    const EBISBox&    a_ebisbox,
				    const IntVectSet& a_cfivs);
};

#endif
