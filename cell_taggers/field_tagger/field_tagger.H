/*!
  @file   field_tagger.H
  @brief  Implementation of cell_tagger that tags cells based on the gradient of the electric field
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _FIELD_TAGGER_
#define _FIELD_TAGGER_

#include "ultralw_tagger.H"

/*!
  @brief Class that tags cells based on the relative curvature of the electric field
  @details This class tets if \f$|\nabla E|\Delta x/|E| > \epsilon_{\textrm{curv}}\f$ OR \f$|E|/E_{\textrm{max}} > \epsilon_{\textrm{mag}}\f$ are fulfilled. If so, the cell is refined. The same tests is used for coarsening, with different epsilons and opposite inequality sign. However, for coarsening, both tests must be fulfilled. 

  Class options
  -------------

      field_tagger.coarsen_curvature = 0.1    # Sets eps_curv for coarsening
      field_tagger.coarsen_magnitude = 0.1    # Sets eps_mag for coarsening 
      field_tagger.refine_curvature  = 0.2    # Sets eps_curv for refinement
      field_tagger.refine_magnitude  = 0.75   # Sets eps_mag for refinement

*/
class field_tagger : public ultralw_tagger {
public:
  
  /*!
    @brief Constructor
  */
  field_tagger();

  /*!
    @brief Destructor
  */
  virtual ~field_tagger();

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
};
#endif
