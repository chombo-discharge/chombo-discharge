/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_CdrPlasmaStreamerTagger.H
  @brief  Implementation of CellTagger that tags cells based on field evaluations or the Townsend coefficient
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaStreamerTagger_H
#define CD_CdrPlasmaStreamerTagger_H

// Our includes
#include <CD_CdrPlasmaFieldTagger.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {
  
    /*!
      @brief Refinement class for morrow_stephens plasma model only
    */
    class CdrPlasmaStreamerTagger : public CdrPlasmaFieldTagger {
    public:
  
      /*!
	@brief Constructor
      */
      CdrPlasmaStreamerTagger();

      /*!
	@brief Constructor
      */
      CdrPlasmaStreamerTagger(const RefCountedPtr<CdrPlasmaPhysics>&     a_physics,
			      const RefCountedPtr<cdr_plasma_stepper>&     a_timeStepper,
			      const RefCountedPtr<AmrMesh>&               a_amr,
			      const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

      /*!
	@brief Destructor
      */
      virtual ~CdrPlasmaStreamerTagger();

      /*!
	@brief Parse class options
      */
      virtual void parseOptions();

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
      virtual bool coarsenCell(const RealVect         a_pos,
			       const Real             a_time,
			       const Real             a_dx,
			       const int              a_lvl,
			       const Vector<Real>     a_tracer,
			       const Vector<RealVect> a_grad_tracer);

      /*!
	@brief Return refinement criterion. 
      */
      virtual bool refineCell(const RealVect         a_pos,
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

#include <CD_NamespaceFooter.H>

#endif
