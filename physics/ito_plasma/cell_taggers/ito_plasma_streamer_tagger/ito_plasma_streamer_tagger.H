/*!
  @file   ito_plasma_streamer_tagger.H
  @brief  Implementation of cell_tagger that tags cells based on the gradient of the electric field
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _ITO_PLASMA_STREAMER_TAGGER_
#define _ITO_PLASMA_STREAMER_TAGGER_

#include "ito_plasma_field_tagger.H"

namespace physics {
  namespace ito_plasma {
  
    /*!
      @brief Refinement class for morrow_stephens plasma model only
    */
    class ito_plasma_streamer_tagger : public ito_plasma_field_tagger {
    public:
  
      /*!
	@brief Constructor
      */
      ito_plasma_streamer_tagger();

      /*!
	@brief Constructor
      */
      ito_plasma_streamer_tagger(const RefCountedPtr<ito_plasma_physics>&     a_physics,
		      const RefCountedPtr<ito_plasma_stepper>&     a_timestepper,
		      const RefCountedPtr<amr_mesh>&               a_amr,
		      const RefCountedPtr<computational_geometry>& a_compgeom);

      /*!
	@brief Destructor
      */
      virtual ~ito_plasma_streamer_tagger();

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

      /*!
	@brief Threshold for coarsening based on curvature
      */
      Real m_coar_curv;

      /*!
	@brief Threshold for refinement based on curvature
      */
      Real m_refi_curv;
  
      /*!
	@brief Threshold for ionization refinement
      */
      Real m_refi_alpha;

      /*!
	@brief Threshold for ionization coarsening
      */
      Real m_coar_alpha;

      /*!
	@brief Set coarsening level
      */
      int m_max_coarsen_level;
    };
  }
}
#endif
