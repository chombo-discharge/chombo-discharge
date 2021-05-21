/*!
  @file   field_tagger.H
  @brief  Implementation of CellTagger that tags cells based on the gradient of the electric field
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _FIELD_TAGGER_
#define _FIELD_TAGGER_

#include "cdr_plasma_field_tagger.H"

#include "CD_NamespaceHeader.H"
namespace Physics {
  namespace CdrPlasma {

    /*!
      @brief Class that tags cells based on the relative curvature of the electric field
    */
    class field_tagger : public cdr_plasma_field_tagger {
    public:
  
      /*!
	@brief Constructor
      */
      field_tagger();

      /*!
	@brief Full constructor
      */
      field_tagger(const RefCountedPtr<CdrPlasmaPhysics>&     a_physics,
		   const RefCountedPtr<cdr_plasma_stepper>&     a_timeStepper,
		   const RefCountedPtr<AmrMesh>&               a_amr,
		   const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

      /*!
	@brief Destructor
      */
      virtual ~field_tagger();

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
  }
}
#include "CD_NamespaceFooter.H"
#endif
