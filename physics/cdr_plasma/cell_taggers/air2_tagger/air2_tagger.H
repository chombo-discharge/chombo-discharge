/*!
  @file   air2_tagger.H
  @brief  Implementation of cell_tagger that tags cells based on the gradient of the electric field
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _AIR2_TAGGER_
#define _AIR2_TAGGER_

#include "ultralw_tagger.H"

/*!
  @brief Refinement class for morrow_zheleznyak plasma model only
*/
class air2_tagger : public ultralw_tagger {
public:
  
  /*!
    @brief Constructor
  */
  air2_tagger();

  /*!
    @brief Destructor
  */
  virtual ~air2_tagger();

  /*!
    @brief Parse class options
  */
  virtual void parse_options();

  /*!
    @brief Return tracer field
  */
  virtual Vector<Real> tracer(const RealVect         a_pos,
			      const Real             a_time,
			      const Real             a_dx,
			      const RealVect         a_E,
			      const Real             a_min_E,
			      const Real             a_max_E,
			      const RealVect         a_grad_E,
			      const Real             a_min_grad_E,
			      const Real             a_max_grad_E);

  /*!
    @brief Return coarsening criterion. 
  */
  virtual bool coarsen_cell(const RealVect         a_pos,
			    const Real             a_time,
			    const Real             a_dx,
			    const int              a_lvl,
			    const Vector<Real>     a_tracer,
			    const Vector<RealVect> a_grad_tracer);

  /*!
    @brief Return refinement criterion. 
  */
  virtual bool refine_cell(const RealVect         a_pos,
			   const Real             a_time,
			   const Real             a_dx,
			   const int              a_lvl,
			   const Vector<Real>     a_tracer,
			   const Vector<RealVect> a_grad_tracer);
protected:

  Real m_coar_curv;
  Real m_refi_curv;
  Real m_refi_alpha;
  Real m_coar_alpha;

  RealVect m_seed_pos;
  Real m_seed_radius;
  Real m_refine_radius;
};
#endif
