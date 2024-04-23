/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
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
{
  CH_TIME("AdvectionDiffusionTagger::AdvectionDiffusionTagger");

  m_solver    = a_solver;
  m_amr       = a_amr;
  m_name      = "AdvectionDiffusion";
  m_verbosity = -1;
  m_realm     = a_solver->getRealm();
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
AdvectionDiffusionTagger::tagCells(EBAMRTags& a_tags)
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
  DataOps::vectorLength(sca, vec);                       // sca = |grad(phi)|
  DataOps::setCoveredValue(sca, 0, 0.0);                 // Set covered cell values to zero.

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

            foundTags = 1;
          }
        }
      };

      // Execute kernel. Regular cells only.
      BoxLoops::loop(cellBox, taggingKernel);
    }
  }

  return (ParallelOps::max(foundTags) > 0) ? true : false;
}

#include <CD_NamespaceFooter.H>
