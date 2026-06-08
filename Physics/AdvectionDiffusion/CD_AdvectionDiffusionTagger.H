/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
   @file   CD_AdvectionDiffusionTagger.H
   @brief  Declaration of the Physics::AdvectionDiffusion::AdvectionDiffusionTagger class
   @author Robert Marskar
*/

#ifndef CD_ADVECTIONDIFFUSIONTAGGER_H
#define CD_ADVECTIONDIFFUSIONTAGGER_H

// Our includes
#include <CD_CellTagger.H>
#include <CD_CdrSolver.H>
#include <CD_AmrMesh.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace AdvectionDiffusion {

    /**
       @brief Cell tagger for the advection-diffusion physics module.
       @details Flags cells for refinement based on two criteria applied to the scalar field phi:
       - Curvature: refine if |grad(phi)|*dx / (phi + epsilon) > m_refCurv.
       - Magnitude: refine only if |phi| > m_refMagn (guards against refining near-zero regions).
       Both thresholds are read from ParmParse at setup time.
    */
    class AdvectionDiffusionTagger : public CellTagger
    {
    public:
      /**
         @brief Disallowed default constructor.
      */
      AdvectionDiffusionTagger() = delete;

      /**
         @brief Full constructor.
         @param[in] a_solver CDR solver whose scalar field is used for tagging decisions.
         @param[in] a_amr    AMR mesh used for gradient computation and level iteration.
      */
      AdvectionDiffusionTagger(RefCountedPtr<CdrSolver>& a_solver, RefCountedPtr<AmrMesh>& a_amr);

      /**
         @brief Destructor.
      */
      virtual ~AdvectionDiffusionTagger();

      /**
         @brief Regrid hook. No operations are performed.
      */
      void
      regrid() override;

      /**
         @brief Parse tagging thresholds from ParmParse.
         @details Reads m_refCurv (curvature threshold) and m_refMagn (magnitude threshold).
      */
      void
      parseOptions() override;

      /**
         @brief Flag cells for refinement based on curvature and magnitude of the scalar field.
         @param[inout] a_tags Cell tags across the AMR hierarchy.
         @return True if at least one cell was tagged on any MPI rank.
      */
      bool
      tagCells(EBAMRTags& a_tags) override;

    protected:
      /**
         @brief Name of the realm on which the solver is registered.
      */
      std::string m_realm;

      /**
         @brief CDR solver providing the scalar field used for tagging.
      */
      RefCountedPtr<CdrSolver> m_solver;

      /**
         @brief AMR mesh used for gradient computation and grid iteration.
      */
      RefCountedPtr<AmrMesh> m_amr;

      /**
         @brief Curvature refinement threshold. Cells with |grad(phi)|*dx/phi above this value are flagged.
      */
      Real m_refCurv;

      /**
         @brief Magnitude refinement threshold. Only cells where |phi| exceeds this value are eligible for flagging.
      */
      Real m_refMagn;
    };
  } // namespace AdvectionDiffusion
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
