/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaStreamerTagger.H
  @brief  Implementation ItoPlasmaTagger which tags on alpha and curvature of field
  @author Robert Marskar
*/

#ifndef CD_ItoPlasmaStreamerTagger_H
#define CD_ItoPlasmaStreamerTagger_H

// Our includes
#include <CD_ItoPlasmaFieldTagger.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoPlasma {

    /*!
      @brief Refinement class for morrow_stephens plasma model only
    */
    template <typename S>
    class ItoPlasmaStreamerTagger : public ItoPlasmaFieldTagger<S>
    {
    public:
      /*!
	@brief Constructor
      */
      ItoPlasmaStreamerTagger();

      /*!
	@brief Constructor
      */
      ItoPlasmaStreamerTagger(const RefCountedPtr<ItoPlasmaPhysics>&      a_physics,
                              const RefCountedPtr<S>&                     a_timeStepper,
                              const RefCountedPtr<AmrMesh>&               a_amr,
                              const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

      /*!
	@brief Destructor
      */
      virtual ~ItoPlasmaStreamerTagger();

      /*!
	@brief Parse class options
      */
      virtual void
      parseOptions();

      /*!
	@brief Parse class options
      */
      virtual void
      parseRuntimeOptions();

      /*!
	@brief Return tracer field
      */
      virtual Vector<Real>
      tracer(const RealVect a_pos,
             const Real     a_time,
             const Real     a_dx,
             const RealVect a_E,
             const Real     a_min_E,
             const Real     a_max_E,
             const RealVect a_grad_E,
             const Real     a_min_grad_E,
             const Real     a_max_grad_E);

      /*!
	@brief Return coarsening criterion. 
      */
      virtual bool
      coarsenCell(const RealVect         a_pos,
                  const Real             a_time,
                  const Real             a_dx,
                  const int              a_lvl,
                  const Vector<Real>     a_tracer,
                  const Vector<RealVect> a_grad_tracer);

      /*!
	@brief Return refinement criterion. 
      */
      virtual bool
      refineCell(const RealVect         a_pos,
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
  } // namespace ItoPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_ItoPlasmaStreamerTaggerImplem.H>

#endif
