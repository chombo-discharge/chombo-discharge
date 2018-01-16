/*!
  @file robinconductivityebbc.H
  @brief Declaration of a Robin type EB boundary condition class for ebconductivityop
  @author Robert Marskar
  @date Jan. 2018
*/

#ifndef _ROBINCONDUCTIVITYEBBC_
#define _ROBINCONDUCTIVITYEBBC_

#include "type_definitions.H"
#include "robin_coef.H"

#include <BaseEBBC.H>

/*!
  @brief Boundary condition class for Robin type boundary conditions on embedded boundaries
*/
class robinconductivityebbc : public ConductivityBaseEBBC {
public:

  /*!
    @brief Constructor
  */
  robinconductivityebbc(const RealVect a_dx, const RealVect a_origin);

  /*!
    @brief Destructor
  */
  ~robinconductivityebbc();

  /*!
    @brief Set bc order
  */
  virtual void set_type(const stencil_type::which_type a_type);

  /*!
    @brief Make it virtual
  */
  virtual void define(const LayoutData<IntVectSet>& a_cfivs, const Real& a_factor);

  /*!
    @brief Set constant coefficients
  */
  virtual void set_coefs(const Real a_aco, const Real b_bco, const Real a_rhs);

  /*!
    @brief Function based coefficients
  */
  virtual void set_coefs(const RefCountedPtr<robin_coef> a_robinco);

  /*!
    @brief Set data-based coefficients
  */
  virtual void set_coefs(const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_aco,
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

  stencil_type::which_type m_type;

  RealVect m_dx;

  RealVect m_origin;

  RefCountedPtr<robin_coef> m_robinco;
  
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

  virtual bool get_taylor_sten(VoFStencil& a_stencil, const VolIndex& a_vof, const EBISBox& a_ebisbox, const IntVectSet& a_cfivs);

  virtual bool get_lsq_sten(VoFStencil& a_stencil, const VolIndex& a_vof, const EBISBox& a_ebisbox, const ProblemDomain& a_domain);
};

#endif
