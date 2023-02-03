/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCStreamerTagger.H
  @brief  Implementation ItoKMCTagger which tags on alpha and curvature of field
  @author Robert Marskar
*/

#ifndef CD_ItoKMCStreamerTagger_H
#define CD_ItoKMCStreamerTagger_H

// Our includes
#include <CD_ItoKMCFieldTagger.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoKMC {

    /*!
      @brief Cell tagging class for streamers.
    */
    template <typename S>
    class ItoKMCStreamerTagger : public ItoKMCFieldTagger<S>
    {
    public:
      /*!
	@brief Constructor. Must subsequently call the define function. 
      */
      ItoKMCStreamerTagger() noexcept;

      /*!
	@brief Full constructor
	@param[in] a_physics Physics implementation. 
	@param[in] a_timeStepper Time integrator. 
	@param[in] a_amr Reference to AmrMesh 
      */
      ItoKMCStreamerTagger(const RefCountedPtr<ItoKMCPhysics>& a_physics,
                           const RefCountedPtr<S>&             a_timeStepper,
                           const RefCountedPtr<AmrMesh>&       a_amr) noexcept;

      /*!
	@brief Destructor
      */
      virtual ~ItoKMCStreamerTagger() noexcept;

      /*!
	@brief Parse class options
      */
      virtual void
      parseOptions() noexcept override;

      /*!
	@brief Parse class options
      */
      virtual void
      parseRuntimeOptions() noexcept override;

    protected:
      /*!
	@brief Threshold for coarsening based on curvature
      */
      Real m_coarsenCurvature;

      /*!
	@brief Threshold for refinement based on curvature
      */
      Real m_refineCurvature;

      /*!
	@brief Threshold for ionization refinement
      */
      Real m_refineAlpha;

      /*!
	@brief Threshold for ionization coarsening
      */
      Real m_coarsenAlpha;

      /*!
	@brief Set coarsening level
      */
      int m_maxCoarsenLevel;

      /*!
	@brief Compute tracer field
	@param[in] a_pos Physical coordinates
	@param[in] a_time Current time
	@param[in] a_dx Grid resolution
	@param[in] a_E Electric field
	@param[in] a_minE Lowest electric field in domain
	@param[in] a_maxE Highest electric field in domain
	@param[in] a_gradE Electric field gradient, i.e. grad(|E|)
	@param[in] a_minGradE Lowest electric field gradient in domain
	@param[in] a_maxGradE Highest electric field gradient in domain
      */
      virtual Vector<Real>
      computeTagFields(const RealVect a_pos,
                       const Real     a_time,
                       const Real     a_dx,
                       const RealVect a_E,
                       const Real     a_minE,
                       const Real     a_maxE,
                       const RealVect a_gradE,
                       const Real     a_minGradE,
                       const Real     a_maxGradE) const noexcept override;

      /*!
	@brief Determine if a particular cell should be coarsened or not. 
	@param[in] a_pos Physical position. 
	@param[in] a_time Current time
	@param[in] a_dx Grid resolution
	@param[in] a_lvl Grid level 
	@param[in] a_tagFields Cell tagging fields
	@param[in] a_gradTagFields Gradient of cell tagging fields
      */
      virtual bool
      coarsenCell(const RealVect         a_pos,
                  const Real             a_time,
                  const Real             a_dx,
                  const int              a_lvl,
                  const Vector<Real>     a_tagFields,
                  const Vector<RealVect> a_gradTagFields) const noexcept override;

      /*!
	@brief Determine if a particular cell should be refined or not. 
	@param[in] a_pos Physical position. 
	@param[in] a_time Current time
	@param[in] a_dx Grid resolution
	@param[in] a_lvl Grid level 
	@param[in] a_tagFields Cell tagging fields
	@param[in] a_gradTagFields Gradient of cell tagging fields
      */
      virtual bool
      refineCell(const RealVect         a_pos,
                 const Real             a_time,
                 const Real             a_dx,
                 const int              a_lvl,
                 const Vector<Real>     a_tagFields,
                 const Vector<RealVect> a_gradTagFields) const noexcept override;
    };
  } // namespace ItoKMC
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_ItoKMCStreamerTaggerImplem.H>

#endif
