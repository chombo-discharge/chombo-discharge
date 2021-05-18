/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file  CD_RobinConductivityEbBc.H
  @brief Declaration of a Robin type EB boundary condition class for EbHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_RobinConductivityEbBc_H
#define CD_RobinConductivityEbBc_H

// Chombo includes
#include <BaseEBBC.H>

// Our includes
#include <CD_IrregStencil.H>
#include <CD_RobinCoefficients.H>
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Boundary condition class for Robin type boundary conditions on embedded boundaries
*/
class RobinConductivityEbBc : public ConductivityBaseEBBC {
public:

  /*!
    @brief Constructor
  */
  RobinConductivityEbBc(const RealVect a_dx, const RealVect a_origin);

  /*!
    @brief Destructor
  */
  ~RobinConductivityEbBc();

  /*!
    @brief Quadrant based stencil or not
  */
  static bool s_quadrant_based;

  /*!
    @brief Radius for least squares
  */
  static int s_lsq_radius;

  /*!
    @brief Set bc order
  */
  virtual void setStencilType(const IrregStencil::StencilType a_type);

  /*!
    @brief Make it virtual
  */
  virtual void define(const LayoutData<IntVectSet>& a_cfivs, const Real& a_factor);

  /*!
    @brief Set constant coefficients
  */
  virtual void setCoefficients(const Real a_aco, const Real b_bco, const Real a_rhs);

  /*!
    @brief Function based coefficients
  */
  virtual void setCoefficients(const RefCountedPtr<RobinCoefficients> a_robinco);

  /*!
    @brief Set data-based coefficients
  */
  virtual void setCoefficients(const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_aco,
			 const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_bco,
			 const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_rhs);

  /*!
    @brief Apply the eb flux
  */
  virtual void applyEBFlux(EBCellFAB&                    a_lphi, 
			   const EBCellFAB&              a_phi, 
			   VoFIterator&                  a_vofit, 
			   const LayoutData<IntVectSet>& a_cfivs, 
			   const DataIndex&              a_dit, 
			   const RealVect&               a_probLo, 
			   const RealVect&               a_dx, 
			   const Real&                   a_factor, 
			   const bool&                   a_useHomogeneous, 
			   const Real&                   a_time);

  /*!
    @brief Get flux stencil
  */
  virtual LayoutData<BaseIVFAB<VoFStencil> >* getFluxStencil(int ivar);

protected:

  LayoutData<BaseIVFAB<VoFStencil> > m_bcstencils;

  IrregStencil::StencilType m_type;

  RealVect m_dx;

  RealVect m_origin;

  RefCountedPtr<RobinCoefficients> m_robinco;
  
  Real m_aco;
  Real m_bco;
  Real m_rhs;

  Real (*m_acofunc)(const RealVect a_pos);
  Real (*m_bcofunc)(const RealVect a_pos);
  Real (*m_rhsfunc)(const RealVect a_pos);

  RefCountedPtr<LevelData<BaseIVFAB<Real> > > m_acodata;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > > m_bcodata;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > > m_rhsdata;

  bool m_const_rhs;
  bool m_data_rhs;
  bool m_func_rhs;

  bool m_const_coeff;
  bool m_data_coeff;
  bool m_func_coeff;

  bool get_taylor_sten(VoFStencil& a_stencil, const VolIndex& a_vof, const EBISBox& a_ebisbox, const IntVectSet& a_cfivs);

  bool get_lsq_sten(VoFStencil& a_stencil, const VolIndex& a_vof, const EBISBox& a_ebisbox, const ProblemDomain& a_domain);
};

#include <CD_NamespaceFooter.H>

#endif
