/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ValidatingBrownianWalkerStepper.H
  @brief  BrownianWalkerStepper subclass that validates the SoA ItoSolver in-situ.
  @author Robert Marskar
  @details
  After each base advance(), deposits the SoA particle weight with each supported deposition
  scheme and prints the max deposited density to pout(). A non-zero, finite max proves the
  SoA deposit pipeline is live on the real cut-cell geometry.
*/

#ifndef CD_VALIDATINGBROWNIANWALKERSTEPPER_H
#define CD_VALIDATINGBROWNIANWALKERSTEPPER_H

// Std includes
#include <algorithm>
#include <cmath>

// Chombo includes
#include <EBCellFAB.H>
#include <FArrayBox.H>
#include <BoxIterator.H>
#include <DataIterator.H>
#include <parstream.H>

// Our includes
#include <CD_BrownianWalkerStepper.H>
#include <CD_ItoSolver.H>
#include <CD_ItoParticle.H>
#include <CD_AmrMesh.H>
#include <CD_ParallelOps.H>
#include <CD_NamespaceHeader.H>

/**
  @brief BrownianWalkerStepper that deposits from the live SoA container each step and prints density stats.
*/
class ValidatingBrownianWalkerStepper : public Physics::BrownianWalker::BrownianWalkerStepper
{
public:
  /**
    @brief Constructor. Forwards the solver to the base stepper.
    @param[in] a_solver Ito solver.
  */
  ValidatingBrownianWalkerStepper(RefCountedPtr<ItoSolver>& a_solver)
    : Physics::BrownianWalker::BrownianWalkerStepper(a_solver)
  {}

  /**
    @brief Advance, then validate the SoA deposition pipeline.
    @param[in] a_dt Time step.
    @return The time step actually taken (from the base stepper).
  */
  Real
  advance(const Real a_dt) override
  {
    const Real ret = Physics::BrownianWalker::BrownianWalkerStepper::advance(a_dt);
    this->validateSoA();
    return ret;
  }

protected:
  /**
    @brief Deposit the live SoA bulk particles and print max deposited density for each scheme.
  */
  void
  validateSoA()
  {
    const Vector<DisjointBoxLayout>& grids       = m_amr->getGrids(m_realm);
    const int                        finestLevel = m_amr->getFinestLevel();

    ParticleContainerSoA<ItoParticle>& soa = m_solver->getParticles(ItoSolver::WhichContainer::Bulk);

    EBAMRCellData mesh;
    m_amr->allocate(mesh, m_realm, m_phase, 1);

    struct Case
    {
      DepositionType       dep;
      CoarseFineDeposition cf;
      const char*          name;
    };
    const Case cases[4] = {{DepositionType::NGP, CoarseFineDeposition::Interp, "NGP/Interp"},
                           {DepositionType::CIC, CoarseFineDeposition::Interp, "CIC/Interp"},
                           {DepositionType::CIC, CoarseFineDeposition::Halo, "CIC/Halo"},
                           {DepositionType::CIC, CoarseFineDeposition::Transition, "CIC/Transition"}};

    for (int c = 0; c < 4; c++) {
      m_amr->depositWeight(mesh, m_realm, m_phase, soa, cases[c].dep, cases[c].cf, false);

      Real maxVal = 0.0;
      for (int lvl = 0; lvl <= finestLevel; lvl++) {
        for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
          const FArrayBox& fab = (*mesh[lvl])[dit()].getFArrayBox();
          for (BoxIterator bit(grids[lvl][dit()]); bit.ok(); ++bit) {
            maxVal = std::max(maxVal, std::abs(fab(bit(), 0)));
          }
        }
      }
      maxVal = ParallelOps::max(maxVal);
      pout() << "  [SoA-validate] " << cases[c].name << "  max deposited density = " << maxVal << std::endl;
    }
  }
};

#include <CD_NamespaceFooter.H>

#endif
