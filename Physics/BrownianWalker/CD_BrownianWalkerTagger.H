/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_BrownianWalkerTagger.H
  @brief  Declaration of the Physics::BrownianWalker::BrownianWalkerTagger class
  @author Robert Marskar
*/

#ifndef CD_BROWNIANWALKERTAGGER_H
#define CD_BROWNIANWALKERTAGGER_H

// Our includes
#include <CD_CellTagger.H>
#include <CD_ItoSolver.H>
#include <CD_AmrMesh.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace BrownianWalker {

    /**
      @brief CellTagger implementation for the Brownian walker physics module.
      @details This tagger does not flag any cells for refinement, but always returns true
      from tagCells so that Driver can trigger load balancing at each regrid step.
    */
    class BrownianWalkerTagger : public CellTagger
    {
    public:
      /**
        @brief Disallowed default constructor.
      */
      BrownianWalkerTagger() = delete;

      /**
        @brief Full constructor.
        @param[in] a_solver Ito solver whose particle data may inform future tagging extensions.
        @param[in] a_amr    AMR mesh providing grid access.
      */
      BrownianWalkerTagger(RefCountedPtr<ItoSolver>& a_solver, RefCountedPtr<AmrMesh>& a_amr);

      /**
        @brief Destructor.
      */
      ~BrownianWalkerTagger();

      /**
        @brief Regrid hook. No operations are performed.
      */
      void
      regrid() override;

      /**
        @brief Parse tagging options from ParmParse. Reads m_refMagn.
      */
      void
      parseOptions() override;

      /**
        @brief Tag cells for refinement. Currently tags no cells but returns true to trigger load balancing.
        @param[inout] a_tags Cell tags across the AMR hierarchy.
        @return Always true.
      */
      bool
      tagCells(EBAMRTags& a_tags) override;

    protected:
      /**
        @brief Ito solver reference, available for future tagging criteria.
      */
      RefCountedPtr<ItoSolver> m_solver;

      /**
        @brief AMR mesh reference.
      */
      RefCountedPtr<AmrMesh> m_amr;

      /**
        @brief Magnitude refinement threshold (reserved for future use).
      */
      Real m_refMagn;
    };
  } // namespace BrownianWalker
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
