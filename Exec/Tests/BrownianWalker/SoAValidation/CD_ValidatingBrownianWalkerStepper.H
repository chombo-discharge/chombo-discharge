/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ValidatingBrownianWalkerStepper.H
  @brief  BrownianWalkerStepper subclass that validates the Dev SoA particle stack in-situ.
  @author Robert Marskar
  @details
  Drops the staged SoA classes (ParticleContainerSoA + EBAMRParticleMeshSoA) in alongside the
  production AoS path of a real BrownianWalker (ItoSolver) run -- real cut-cell geometry
  (RodDielectric) and real AMR. After each base advance(), it mirrors the live bulk ItoParticles
  (position + weight) into a SoA container on the SAME realm/grids, deposits the weight with BOTH the
  production EBAMRParticleMesh and the SoA EBAMRParticleMeshSoA, and prints max|AoS-SoA| per
  deposition scheme. This is a validation harness only -- it lives in the app, not the library, so the
  production module is untouched.
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
#include <List.H>
#include <RealVect.H>
#include <parstream.H>

// Our includes
#include <CD_BrownianWalkerStepper.H>
#include <CD_ItoSolver.H>
#include <CD_ItoParticle.H>
#include <CD_AmrMesh.H>
#include <CD_EBAMRParticleMesh.H>
#include <CD_ParallelOps.H>
#include <CD_ParticleContainerSoA.H>
#include <CD_EBAMRParticleMeshSoA.H>
#include <CD_NamespaceHeader.H>

/**
  @brief BrownianWalkerStepper that additionally runs an AoS-vs-SoA deposition parity check each step.
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
    @brief Advance, then validate the SoA deposition against production.
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
    @brief Mirror the live bulk particles into a SoA container and compare AoS vs SoA deposition.
  */
  void
  validateSoA()
  {
    const Vector<DisjointBoxLayout>&          grids          = m_amr->getGrids(m_realm);
    const Vector<ProblemDomain>&              domains        = m_amr->getDomains();
    const Vector<Real>&                       dx             = m_amr->getDx();
    const Vector<int>&                        refRat         = m_amr->getRefinementRatios();
    const RealVect                            probLo         = m_amr->getProbLo();
    const int                                 finestLevel    = m_amr->getFinestLevel();
    const int                                 blockingFactor = m_amr->getBlockingFactor();
    const int                                 nghost         = m_amr->getNumberOfGhostCells();
    const Vector<RefCountedPtr<EBLevelGrid>>& eblgs          = m_amr->getEBLevelGrid(m_realm, m_phase);

    // Live production particles (already remapped + deposited this step).
    ParticleContainer<ItoParticle>& aos = m_solver->getParticles(ItoSolver::WhichContainer::Bulk);

    // Mirror (position, weight) into a SoA container on the SAME grids (identical placement).
    ParticleContainerSoA<> soa;
    soa.define(grids, domains, dx, refRat, probLo, blockingFactor, finestLevel, m_realm);
    for (int lvl = 0; lvl <= finestLevel; lvl++) {
      for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
        const List<ItoParticle>&      L    = aos[lvl][dit()].listItems();
        ParticleContainerSoA<>::Leaf& leaf = soa[lvl][dit()];
        for (ListIterator<ItoParticle> lit(L); lit.ok(); ++lit) {
          leaf.append(lit().position(), lit().weight(), NoPayload{});
        }
      }
    }

    EBAMRParticleMeshSoA pmSoA;
    pmSoA.define(eblgs, refRat, dx, probLo, nghost, finestLevel);
    EBAMRParticleMesh& pmAoS = m_amr->getParticleMesh(m_realm, m_phase);

    EBAMRCellData meshAoS;
    EBAMRCellData meshSoA;
    m_amr->allocate(meshAoS, m_realm, m_phase, 1);
    m_amr->allocate(meshSoA, m_realm, m_phase, 1);

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
      pmAoS.deposit<ItoParticle, const Real&, &ItoParticle::weight>(meshAoS, aos, cases[c].dep, cases[c].cf, false);
      pmSoA.depositWeight(meshSoA, soa, cases[c].dep, cases[c].cf, false);

      Real maxAbs = 0.0;
      for (int lvl = 0; lvl <= finestLevel; lvl++) {
        for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
          const FArrayBox& a = (*meshAoS[lvl])[dit()].getFArrayBox();
          const FArrayBox& b = (*meshSoA[lvl])[dit()].getFArrayBox();
          for (BoxIterator bit(grids[lvl][dit()]); bit.ok(); ++bit) {
            maxAbs = std::max(maxAbs, std::abs(a(bit(), 0) - b(bit(), 0)));
          }
        }
      }
      maxAbs = ParallelOps::max(maxAbs);
      pout() << "  [SoA-validate] " << cases[c].name << "  max|AoS-SoA| = " << maxAbs << endl;
    }
  }
};

#include <CD_NamespaceFooter.H>

#endif
