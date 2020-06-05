/*!
  @file   ito_plasma_field_tagger.H
  @brief  Declaration of an abstract field-only tagging class for ito plasmas
  @author Robert Marskar
  @date   May. 2018
*/

#ifndef _ITO_PLASMA_FIELD_TAGGER_
#define _ITO_PLASMA_FIELD_TAGGER_

#include "ito_plasma_tagger.H"

namespace physics {
  namespace ito_plasma {
  
    /*!
      @brief Lightweight tagging class. Only tags cells based on electric field. 
      @details This class fetches the point-wise electric field. The class is pure, so users must implement 
      the actual tagging criterion themselves. 
    */
    class ito_plasma_field_tagger : public ito_plasma_tagger {
    public:

      /*!
	@brief Weak constructor, must subsequently call define
      */
      ito_plasma_field_tagger();

      /*!
	@brief Destructor
      */
      ~ito_plasma_field_tagger();

      /*!
	@brief Compute tracer fields
      */
      void compute_tracers();

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
	@brief Storage for electric field
      */
      EBAMRCellData m_E;

      /*!
	@brief Storage for electric field
      */
      EBAMRCellData m_grad_E;

      /*!
	@brief Allocate memory
      */
      void allocate_storage();

      /*!
	@brief Deallocate memory
      */
      void deallocate_storage();

      /*!
	@brief Compute electric field
      */
      void compute_E(EBAMRCellData& a_E, EBAMRCellData& a_grad_E);
    
      /*!
	@brief Compute tracer field
      */
      virtual Vector<Real> tracer(const RealVect          a_pos,
				  const Real              a_time,
				  const Real              a_dx,
				  const RealVect          a_E,
				  const Real              a_min_E,
				  const Real              a_max_E,
				  const RealVect          a_grad_E,
				  const Real              a_min_grad_E,
				  const Real              a_max_grad_E) = 0;


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
  }
}
#endif
