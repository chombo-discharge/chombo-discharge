/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ItoKMCStreamerTagger.H
  @brief  Declaration of the Physics::ItoKMC::ItoKMCStreamerTagger CellTagger
  @author Robert Marskar
*/

#ifndef CD_ITOKMCSTREAMERTAGGER_H
#define CD_ITOKMCSTREAMERTAGGER_H

// Our includes
#include <CD_ItoKMCFieldTagger.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoKMC {

    /**
      @brief CellTagger for the ItoKMC module that refines on the ionization coefficient and field curvature.
      @details Refinement is triggered when @f$ \alpha \cdot \Delta x @f$ exceeds a user-specified threshold or
      when the normalised curvature of |E| exceeds another threshold. Coarsening uses the same fields with
      separate (lower) thresholds.
    */
    template <typename S>
    class ItoKMCStreamerTagger : public ItoKMCFieldTagger<S>
    {
    public:
      /**
	@brief Constructor. Must subsequently call the define function. 
      */
      ItoKMCStreamerTagger() noexcept;

      /**
	@brief Full constructor
	@param[in] a_physics Physics implementation. 
	@param[in] a_timeStepper Time integrator. 
	@param[in] a_amr Reference to AmrMesh 
      */
      ItoKMCStreamerTagger(const RefCountedPtr<ItoKMCPhysics>& a_physics,
                           const RefCountedPtr<S>&             a_timeStepper,
                           const RefCountedPtr<AmrMesh>&       a_amr) noexcept;

      /**
	@brief Destructor
      */
      virtual ~ItoKMCStreamerTagger() noexcept;

      /**
	@brief Parse class options
      */
      virtual void
      parseOptions() noexcept override;

      /**
	@brief Parse runtime-configurable class options.
      */
      virtual void
      parseRuntimeOptions() noexcept override;

    protected:
      /**
	@brief Threshold for coarsening based on curvature
      */
      Real m_coarsenCurvature;

      /**
	@brief Threshold for refinement based on curvature
      */
      Real m_refineCurvature;

      /**
	@brief Threshold for ionization refinement
      */
      Real m_refineAlpha;

      /**
	@brief Threshold for ionization coarsening
      */
      Real m_coarsenAlpha;

      /**
	@brief Set coarsening level
      */
      int m_maxCoarsenLevel;

      /**
	@brief Compute per-cell tagging fields from the electric field and its gradient.
	@param[in] a_pos      Physical coordinates.
	@param[in] a_time     Current simulation time.
	@param[in] a_dx       Grid resolution.
	@param[in] a_E        Electric field vector at this cell.
	@param[in] a_minE     Minimum |E| over the domain.
	@param[in] a_maxE     Maximum |E| over the domain.
	@param[in] a_gradE    grad(|E|) at this cell.
	@param[in] a_minGradE Minimum |grad(|E|)| over the domain.
	@param[in] a_maxGradE Maximum |grad(|E|)| over the domain.
	@return Per-cell tagging fields [normalised curvature, normalised alpha*dx].
      */
      virtual Vector<Real>
      computeTagFields(const RealVect& a_pos,
                       const Real      a_time,
                       const Real      a_dx,
                       const RealVect& a_E,
                       const Real      a_minE,
                       const Real      a_maxE,
                       const RealVect& a_gradE,
                       const Real      a_minGradE,
                       const Real      a_maxGradE) const noexcept override;

      /**
	@brief Determine if a particular cell should be coarsened.
	@param[in] a_pos           Physical position.
	@param[in] a_time          Current simulation time.
	@param[in] a_dx            Grid resolution.
	@param[in] a_lvl           Grid level.
	@param[in] a_tagFields     Cell tagging field values.
	@param[in] a_gradTagFields Gradient of cell tagging fields.
	@return True if the cell should be coarsened.
      */
      virtual bool
      coarsenCell(const RealVect&         a_pos,
                  const Real              a_time,
                  const Real              a_dx,
                  const int               a_lvl,
                  const Vector<Real>&     a_tagFields,
                  const Vector<RealVect>& a_gradTagFields) const noexcept override;

      /**
	@brief Determine if a particular cell should be refined.
	@param[in] a_pos           Physical position.
	@param[in] a_time          Current simulation time.
	@param[in] a_dx            Grid resolution.
	@param[in] a_lvl           Grid level.
	@param[in] a_tagFields     Cell tagging field values.
	@param[in] a_gradTagFields Gradient of cell tagging fields.
	@return True if the cell should be refined.
      */
      virtual bool
      refineCell(const RealVect&         a_pos,
                 const Real              a_time,
                 const Real              a_dx,
                 const int               a_lvl,
                 const Vector<Real>&     a_tagFields,
                 const Vector<RealVect>& a_gradTagFields) const noexcept override;
    };
  } // namespace ItoKMC
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_ItoKMCStreamerTaggerImplem.H>

#endif
