/*!
  @file   gradient_tagger.H
  @brief  Implementation of cell_tagger that tags cells based on the gradient of the electric field
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _GRADIENT_TAGGER_
#define _GRADIENT_TAGGER_

#include "full_tagger.H"

/*!
  @brief Class that tags cells based on the relative curvature of the cdr equation densities
*/
class gradient_tagger : public full_tagger {
public:
  
  /*!
    @brief Constructor
  */
  gradient_tagger();

  /*!
    @brief Destructor
  */
  virtual ~gradient_tagger();

  /*!
    @brief Define function. 
    @details I'm overriding this so I can fetch the number of species from plasma kinetics
  */
  virtual void define(const RefCountedPtr<plasma_kinetics>&        a_plaskin,
		      const RefCountedPtr<time_stepper>&           a_timestepper,
		      const RefCountedPtr<AmrMesh>&               a_amr,
		      const RefCountedPtr<computational_geometry>& a_compgeom,
		      const RefCountedPtr<physical_domain>&        a_physdom);

  /*!
    @brief Parse options
  */
  virtual void parseOptions();

  /*!
    @brief Build tracer field
  */
  virtual Vector<Real> tracer(const RealVect         a_pos,
			      const Real             a_time,
			      const Real             a_dx,
			      const RealVect         a_E,
			      const Real             a_min_E,
			      const Real             a_max_E,
			      const RealVect         a_grad_E,
			      const Real             a_min_grad_E,
			      const Real             a_max_grad_E,
			      const Real             a_rho,
			      const Real             a_min_rho,
			      const Real             a_max_rho,
			      const RealVect         a_grad_rho,
			      const Real             a_min_grad_rho,
			      const Real             a_max_grad_rho,
			      const Vector<Real>     a_ion_densities,
			      const Vector<Real>     a_min_ion_densities,
			      const Vector<Real>     a_max_ion_densities,
			      const Vector<RealVect> a_ion_gradients,
			      const Vector<Real>     a_min_ion_gradients,
			      const Vector<Real>     a_max_ion_gradients,
			      const Vector<Real>     a_photon_densities,
			      const Vector<Real>     a_min_photon_densities,
			      const Vector<Real>     a_max_photon_densities);

  /*!
    @brief Coarsen a cell based on a tracer field
  */
  virtual bool coarsen_cell(const RealVect         a_pos,
			    const Real             a_time,
			    const Real             a_dx,
			    const int              a_lvl,
			    const Vector<Real>     a_tracer,
			    const Vector<RealVect> a_grad_tracer);

  /*!
    @brief Refine a cell based on a tracer field
  */
  virtual bool refine_cell(const RealVect         a_pos,
			   const Real             a_time,
			   const Real             a_dx,
			   const int              a_lvl,
			   const Vector<Real>     a_tracer,
			   const Vector<RealVect> a_grad_tracer);

protected:

  /*!
    @brief Threshold for coarsening based on curvature
  */
  Real m_coar_curv;

  /*!
    @brief Threshold for refinement based on curvature
  */
  Real m_refi_curv;

  /*!
    @brief Relative cut-off. If magnitude falls below this, don't refine. 
  */
  Real m_cutoff;
#include "CD_NamespaceFooter.H"
#endif
