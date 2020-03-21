/*!
  @file   full_tagger.H
  @brief  Class for tagging cells across an AMR hierarchy
  @author Robert Marskar
  @date   Nov. 2017
  @todo   Rename to full_tagger
*/

#ifndef _FULL_TAGGER_
#define _FULL_TAGGER_

#include "cdr_plasma_tagger.H"

/*!
  @brief Abstract class for tagging cells across an AMR hierarchy
*/
class full_tagger : public cdr_plasma_tagger {
public:

  /*!
    @brief Weak constructor, must subsequently call define
  */
  full_tagger();

  /*!
    @brief Destructor
  */
  virtual ~full_tagger();

  /*!
    @brief Compute tracer fields
  */
  virtual void compute_tracers();

  /*!
    @brief Parse options
  */
  virtual void parse_options() = 0;
  
protected:

  /*!
    @brief Scratch storage
  */
  EBAMRCellData m_scratch;

  /*!
    @brief Electric field
  */
  EBAMRCellData m_E;

  /*!
    @brief grad(|E|)
  */
  EBAMRCellData m_grad_E;

  /*!
    @brief Space charge density
  */
  EBAMRCellData m_rho;

  /*!
    @brief grad(rho)
  */
  EBAMRCellData m_grad_rho;

  /*!
    @brief
  */
  Vector<EBAMRCellData> m_cdr_densities;

  /*!
    @brief
  */
  Vector<EBAMRCellData> m_cdr_gradients;

  /*!
    @brief
  */
  Vector<EBAMRCellData> m_rte_densities;

  /*!
    @brief Allocate storage
  */
  virtual void allocate_storage();

  /*!
    @brief Deallocate storage. 
  */
  virtual void deallocate_storage();

  /*!
    @brief Compute centroid-centered ion densities
  */
  virtual void compute_cdr_densities(Vector<EBAMRCellData>& a_cdr_densities);

  /*!
    @brief Compute centroid-centered ion densities
  */
  virtual void compute_cdr_gradients(Vector<EBAMRCellData>& a_cdr_gradients);

  /*!
    @brief Compute electric field and gradient of |E| at centroids
  */
  virtual void compute_E(EBAMRCellData& a_E, EBAMRCellData& a_grad_E);

  /*!
    @brief Compute space charge and its gardient
  */
  virtual void compute_rho(EBAMRCellData& a_rho, EBAMRCellData& a_grad_rho);

  /*!
    @brief Compute the rte densities
  */
  virtual void compute_rte_densities(Vector<EBAMRCellData>& a_rte_densities);

  /*!
    @brief Build tracer field
    @param[in] a_pos Position
    @param[in] a_time Time
    @param[in] a_dx Level resolution
    @param[in] a_E Electric field
    @param[in] a_min_E Minimum electric field in domain
    @param[in] a_max_E Maximum electric field in domain
    @param[in] a_grad_E grad(|E|)
    @param[in] a_max_grad_E Maximum grad(|E|)
    @param[in] a_min_grad_E Maximum grad(|E|)
    @param[in] a_rho Space charge density
    @param[in] a_min_rho Minimum rho
    @param[in] a_max_rho Maximum rho
    @param[in] a_grad_rho grad(rho)
    @param[in] a_max_grad_rho Maximum grad(rho)
    @param[in] a_min_grad_rho Minimum  grad(rho)
    @param[in] a_ion_densities Ion densities
    @param[in] a_min_ion_densities Lowest densities in domain
    @param[in] a_max_ion_densities Highest
    @param[in] a_ion_gradients Ion gradients
    @param[in] a_min_ion_gradients Smallest ion gradients in domain
    @param[in] a_max_ion_gradients Largest ion gradients in domain
    @param[in] a_photon_densities Photon densities
    @param[in] a_min_photon_densities Minium photon densities in domain
    @param[in] a_max_photon_densities Maximum photon densities in domain
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
			      const Vector<Real>     a_max_photon_densities) = 0;

  /*!
    @brief Coarsen a cell based on a tracer field
  */
  virtual bool coarsen_cell(const RealVect         a_pos,
			    const Real             a_time,
			    const Real             a_dx,
			    const int              a_lvl,
			    const Vector<Real>     a_tracer,
			    const Vector<RealVect> a_grad_tracer) = 0;

  /*!
    @brief Refine a cell based on a tracer field
  */
  virtual bool refine_cell(const RealVect         a_pos,
			   const Real             a_time,
			   const Real             a_dx,
			   const int              a_lvl,
			   const Vector<Real>     a_tracer,
			   const Vector<RealVect> a_grad_tracer) = 0;
};

#endif
