/*!
  @file   advection_kinetics2.H
  @brief  Declaration of the Morrow-Lowke plasma model for air
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _ADVECTION_KINETICS2_
#define _ADVECTION_KINETICS2_

#include "plasma_kinetics.H"
#include "perlin_if.H"

/*!
  @brief   Class the implements advection-diffusion of a single species
*/
class advection_kinetics2 : public plasma_kinetics {
public:
  
  advection_kinetics2();

  virtual ~advection_kinetics2();
  
  virtual Vector<Real> compute_cdr_diffusion_coefficients(const Real&         a_time,
							  const RealVect&     a_pos,
							  const RealVect&     a_E,
							  const Vector<Real>& a_cdr_densities) const;

  virtual Vector<RealVect> compute_cdr_velocities(const Real&         a_time,
						  const RealVect&     a_pos,
						  const RealVect&     a_E,
						  const Vector<Real>& a_cdr_densities) const;

  virtual Vector<Real> compute_cdr_source_terms(const Real              a_time,
						const RealVect&         a_pos,
						const RealVect&         a_E,
						const RealVect&         a_gradE,
						const Vector<Real>&     a_cdr_densities,
						const Vector<Real>&     a_rte_densities,
						const Vector<RealVect>& a_grad_cdr) const;

  virtual Vector<Real> compute_cdr_electrode_fluxes(const Real&         a_time,
						    const RealVect&     a_pos,
						    const RealVect&     a_normal,
						    const RealVect&     a_E,
						    const Vector<Real>& a_cdr_densities,
						    const Vector<Real>& a_cdr_velocities,
						    const Vector<Real>& a_cdr_gradients,
						    const Vector<Real>& a_rte_fluxes,
						    const Vector<Real>& a_extrap_cdr_fluxes) const;

  virtual Vector<Real> compute_cdr_dielectric_fluxes(const Real&         a_time,
						     const RealVect&     a_pos,
						     const RealVect&     a_normal,
						     const RealVect&     a_E,
						     const Vector<Real>& a_cdr_densities,
						     const Vector<Real>& a_cdr_velocities,
						     const Vector<Real>& a_cdr_gradients,
						     const Vector<Real>& a_rte_fluxes,
						     const Vector<Real>& a_extrap_cdr_fluxes) const;

  virtual Vector<Real> compute_rte_source_terms(const Real&         a_time,
						const RealVect&     a_pos,
						const RealVect&     a_E,
						const Vector<Real>& a_cdr_densities) const;

  virtual Real initial_sigma(const Real a_time, const RealVect& a_pos) const;

  /*!
    @brief Advected class
  */
  class phi_advect : public species {
  public:

    phi_advect();

    ~phi_advect();

    Real initial_data(const RealVect a_pos, const Real a_time) const;
    
  protected:

    RealVect m_center;
    RealVect m_width;
  };

  /*!
    @brief The other advected class
  */
  class psi_advect : public species {
  public:

    psi_advect();

    ~psi_advect();

    Real initial_data(const RealVect a_pos, const Real a_time) const;
  };

protected:

  bool m_outflow;

  Real m_beta;
};

#endif
