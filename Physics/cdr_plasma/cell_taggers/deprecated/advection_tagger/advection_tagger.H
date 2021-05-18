/*!
  @file   gradient_tagger.H
  @brief  Implementation of CellTagger that tags cells based on the gradient of the electric field
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _GRADIENT_TAGGER_
#define _GRADIENT_TAGGER_

#include "ultralw_tagger.H"

/*!
  @brief Class that tags cells based on the relative curvature of the electric field
  @details This class tets if \f$|\nabla E|\Delta x/|E| > \epsilon_{\textrm{curv}}\f$ OR \f$|E|/E_{\textrm{max}} > \epsilon_{\textrm{mag}}\f$ are fulfilled. If so, the cell is refined. The same tests is used for coarsening, with different epsilons and opposite inequality sign. However, for coarsening, both tests must be fulfilled. 

  Class options
  -------------

  gradient_tagger.coarsen_curvature = 0.1    # Sets eps_curv for coarsening
  gradient_tagger.coarsen_magnitude = 0.1    # Sets eps_mag for coarsening 
  gradient_tagger.refine_curvature  = 0.2    # Sets eps_curv for refinement
  gradient_tagger.refine_magnitude  = 0.75   # Sets eps_mag for refinement

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
    @brief Return tracer field
  */
  virtual Vector<Real> tracer(const RealVect&         a_pos,
			      const Real&             a_time,
			      const Real&             a_dx,
			      const RealVect&         a_E,
			      const Real&             a_min_E,
			      const Real&             a_max_E,
			      const RealVect&         a_grad_E,
			      const Real&             a_min_grad_E,
			      const Real&             a_max_grad_E);

  /*!
    @brief Return coarsening criterion. 
  */
  virtual bool coarsen_cell(const RealVect&         a_pos,
			    const Real&             a_time,
			    const Real&             a_dx,
			    const int&              a_lvl,
			    const Vector<Real>&     a_tracer,
			    const Vector<RealVect>& a_grad_tracer);

  /*!
    @brief Return refinement criterion. 
  */
  virtual bool refine_cell(const RealVect&         a_pos,
			   const Real&             a_time,
			   const Real&             a_dx,
			   const int&              a_lvl,
			   const Vector<Real>&     a_tracer,
			   const Vector<RealVect>& a_grad_tracer);
protected:

  /*!
    @brief Threshold for coarsening based on curvature
  */
  Real m_coar_curv;

  /*!
    @brief Threshold for coarsening based on magnitude
  */
  Real m_coar_mag;

  /*!
    @brief Threshold for refinement based on curvature
  */
  Real m_refi_curv;
  
  /*!
    @brief Threshold for coarsening based on magnitude
  */
  Real m_refi_mag;
#include "CD_NamespaceFooter.H"
#endif
