/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
   @file   CD_AdvectionDiffusionTagger.cpp
   @brief  Implementation of CD_AdvectionDiffusionTagger.H
   @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_AdvectionDiffusionTagger.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::AdvectionDiffusion;

AdvectionDiffusionTagger::AdvectionDiffusionTagger(RefCountedPtr<CdrSolver>& a_solver, RefCountedPtr<AmrMesh>& a_amr)
  : m_realm(a_solver->getRealm()), m_solver(a_solver), m_amr(a_amr)
{
  CH_TIME("AdvectionDiffusionTagger::AdvectionDiffusionTagger");

  m_name      = "AdvectionDiffusion";
  m_verbosity = -1;
}

AdvectionDiffusionTagger::~AdvectionDiffusionTagger()
{
  CH_TIME("AdvectionDiffusionTagger::~AdvectionDiffusionTagger");
}

void
AdvectionDiffusionTagger::regrid()
{
  CH_TIME("AdvectionDiffusionTagger::regrid");
}

void
AdvectionDiffusionTagger::parseOptions()
{
  CH_TIME("AdvectionDiffusionTagger::parseOptions");

  ParmParse pp(m_name.c_str());

  pp.get("refine_curv", m_refCurv);
  pp.get("refine_magn", m_refMagn);

  this->parseBuffer(); // Derived from CellTagger -- sets buffer region for cell refinement.
}

bool
AdvectionDiffusionTagger::tagCells(EBAMRTags& a_tags) // NOLINT(readability-convert-member-functions-to-static)
{
  CH_TIME("AdvectionDiffusionTagger::tagCells");

  // TLDR: We allocate storage for computing grad(phi) and |grad(phi)|. We then run through the cells on each level and
  //       refine cells where |grad(phi)|*dx/phi > m_refCurv or where phi > m_refMagn.

  // Data allocation
  EBAMRCellData sca;
  EBAMRCellData vec;

  m_amr->allocate(sca, m_realm, phase::gas, 1);
  m_amr->allocate(vec, m_realm, phase::gas, SpaceDim);

  // Get cell-centered state from the solver.
  const EBAMRCellData& state = m_solver->getPhi();

  // Update ghost cells so we can safely compute gradient.
  DataOps::copy(sca, state);
  m_amr->interpGhostMG(sca, m_realm, phase::gas);

  // Compute the gradient
  m_amr->computeGradient(vec, sca, m_realm, phase::gas); // vec =  grad(phi)
  DataOps::vectorLength(sca,
                        vec,
                        m_amr->getNotCoveredCells(m_realm, phase::gas),
                        m_amr->getMultiCutVofIterator(m_realm, phase::gas));
  DataOps::setCoveredValue(sca,
                           m_amr->getCoveredCells(m_realm, phase::gas),
                           0,
                           0.0); // Set covered cell values to zero.

  int foundTags = 0;

  // Never tag on finest possible AMR level.
  const int finestLevel    = m_amr->getFinestLevel();
  const int maxLevel       = m_amr->getMaxAmrDepth();
  const int finestTagLevel = (finestLevel == maxLevel) ? maxLevel - 1 : finestLevel;

  for (int lvl = 0; lvl <= finestTagLevel; lvl++) {
    DataOps::scale(*sca[lvl], m_amr->getDx()[lvl]); // sca = |grad(phi)|*dx

    constexpr Real SAFETY = 1.E-6; // To avoid division by zero

    // Get this on this level
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(+ : foundTags)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box        cellBox = dbl[din];
      const EBCellFAB& gradPhi = (*sca[lvl])[din];   // |grad(phi)|
      const EBCellFAB& phi     = (*state[lvl])[din]; //  phi
      const EBISBox&   ebisBox = phi.getEBISBox();

      // Clear all previous cell tags.
      DenseIntVectSet& tags = (*a_tags[lvl])[din];
      tags.makeEmptyBits();

      // Go through regular grid cells.
      const BaseFab<Real>& gradReg = gradPhi.getSingleValuedFAB();
      const BaseFab<Real>& phiReg  = phi.getSingleValuedFAB();

      // Kernel for tagging
      auto taggingKernel = [&](const IntVect& iv) -> void {
        if (ebisBox.isRegular(iv)) {

          const Real curv = std::abs(gradReg(iv, 0)) / (SAFETY + std::abs(phiReg(iv, 0)));
          if (curv > m_refCurv && std::abs(phiReg(iv, 0)) > m_refMagn) {
            tags |= iv;

            foundTags += 1;
          }
        }
      };

      // Execute kernel. Regular cells only. Not vectorizable: data-dependent curvature test + DenseIntVectSet
      // insertion (tags |= iv). reduction(+ : foundTags) + per-box tags -> no race. One-time tagging.
      BoxLoops::loop<D_DECL(1, 1, 1)>(cellBox, taggingKernel);
    }
  }

  return (ParallelOps::max(foundTags) > 0) ? true : false;
}

#include <CD_NamespaceFooter.H>
