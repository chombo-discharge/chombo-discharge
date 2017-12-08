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
    @brief Destructor
  */
  virtual ~jump_bc();

  /*!
    @brief Define function
  */
  virtual void define(const MFLevelGrid& a_mflg, const Real& a_dx, const int a_order);

  /*!
    @brief Build stencils for matching
  */
  virtual void build_stencils();

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
    @brief Get a second order stencil for the current vof
  */
  virtual bool get_second_order_sten(Real& a_weight, VoFStencil& a_stencil, const VolIndex& a_vof, const EBISBox& a_ebisbox);

  /*!
    @brief Get a first order stencil
  */
  virtual bool get_first_order_sten(Real& a_weight, VoFStencil& a_stencil, const VolIndex& a_vof, const EBISBox& a_ebisbox);
};

#endif
