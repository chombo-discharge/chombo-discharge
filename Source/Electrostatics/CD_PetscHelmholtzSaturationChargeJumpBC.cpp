/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzSaturationChargeJumpBC.cpp
 * @brief  Implementation of CD_PetscHelmholtzSaturationChargeJumpBC.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Std includes
#include <algorithm>
#include <vector>

// Chombo includes
#include <CH_Timer.H>
#include <BaseIVFactory.H>

// Our includes
#include <CD_PetscHelmholtzSaturationChargeJumpBC.H>
#include <CD_NamespaceHeader.H>

PetscHelmholtzSaturationChargeJumpBC::PetscHelmholtzSaturationChargeJumpBC(const phase::which_phase a_phase) noexcept
  : m_phase(a_phase)
{
  CH_TIME("PetscHelmholtzSaturationChargeJumpBC::PetscHelmholtzSaturationChargeJumpBC");
}

PetscHelmholtzSaturationChargeJumpBC::~PetscHelmholtzSaturationChargeJumpBC() noexcept
{
  CH_TIME("PetscHelmholtzSaturationChargeJumpBC::~PetscHelmholtzSaturationChargeJumpBC");
  if (m_verbose) {
    pout() << "PetscHelmholtzSaturationChargeJumpBC::~PetscHelmholtzSaturationChargeJumpBC" << endl;
  }
}

void
PetscHelmholtzSaturationChargeJumpBC::clear() noexcept
{
  CH_TIME("PetscHelmholtzSaturationChargeJumpBC::clear");
  if (m_verbose) {
    pout() << "PetscHelmholtzSaturationChargeJumpBC::clear" << endl;
  }

  m_sigma.clear();

  PetscHelmholtzJumpBC::clear();
}

void
PetscHelmholtzSaturationChargeJumpBC::getInterfaceFluxStencil(Real&            a_dataWeight,
                                                              PetscStencil&    a_stencil,
                                                              const int        a_phase,
                                                              const int        a_level,
                                                              const DataIndex& a_din,
                                                              const IntVect&   a_iv,
                                                              const Real       a_bCoefPhase0,
                                                              const Real       a_bCoefPhase1) const noexcept
{
  CH_TIME("PetscHelmholtzSaturationChargeJumpBC::getInterfaceFluxStencil");
  if (m_verbose) {
    pout() << "PetscHelmholtzSaturationChargeJumpBC::getInterfaceFluxStencil" << endl;
  }

  CH_assert(m_isDefined);
  CH_assert(a_phase == 0 || a_phase == 1);

  a_dataWeight = 0.0;
  a_stencil.clear();

  // No flux into the side of m_phase.
  if (a_phase == m_phase) {
    return;
  }

  // D_z = w_z * phi_B + S_z = 0 on the side z = m_phase gives phi_B = -S_z/w_z, so that on the other side
  //
  //   b_o D_o = b_o (w_o phi_B + S_o) = b_o S_o - b_o w_o / w_z S_z.
  const int z = m_phase;
  const int o = a_phase;

  Real         w[2];
  PetscStencil S[2];

  for (int k = 0; k < 2; k++) {
    this->getNormalDerivativeStencil(w[k], S[k], k, a_level, a_din, VolIndex(a_iv, 0));
  }

  // Without a stencil on side z, the normal derivative is zero on the other side too, as in FieldSolverGMG: there is no
  // flux through the interface.
  if (w[z] == 0.0) {
    return;
  }

  const Real b[2] = {a_bCoefPhase0, a_bCoefPhase1};

  PetscStencil stencilThis  = S[o];
  PetscStencil stencilOther = S[z];

  stencilThis *= b[o];
  stencilOther *= -b[o] * w[o] / w[z];

  a_stencil = stencilThis;
  a_stencil += stencilOther;
  a_stencil.compress();
}

Real
PetscHelmholtzSaturationChargeJumpBC::getJumpData(const int /*a_level*/,
                                                  const DataIndex& /*a_din*/,
                                                  const IntVect& /*a_iv*/,
                                                  const EBAMRIVData& /*a_sigma*/) const noexcept
{
  CH_TIME("PetscHelmholtzSaturationChargeJumpBC::getJumpData");
  if (m_verbose) {
    pout() << "PetscHelmholtzSaturationChargeJumpBC::getJumpData" << endl;
  }

  return 0.0;
}

void
PetscHelmholtzSaturationChargeJumpBC::updateJumpData(const Vec& a_phi, const MFAMRIVData& a_bCoefIrreg) noexcept
{
  CH_TIME("PetscHelmholtzSaturationChargeJumpBC::updateJumpData");
  if (m_verbose) {
    pout() << "PetscHelmholtzSaturationChargeJumpBC::updateJumpData" << endl;
  }

  CH_assert(m_isDefined);

  const Vector<RefCountedPtr<MFLevelGrid>>& mflg = m_petscGrid->getMFLevelGrids();

  const int  numLevels    = mflg.size();
  const bool hasTwoPhases = (m_petscGrid->getNumPhases() == 2);
  const int  z            = m_phase;
  const int  o            = 1 - m_phase;

  // sigma lives on the cut cells of the first phase, and is zero except on the interface cells.
  m_sigma.resize(numLevels);

  for (int lvl = 0; lvl < numLevels; lvl++) {
    const EBLevelGrid&       eblg  = mflg[lvl]->getEBLevelGrid(0);
    const EBISLayout&        ebisl = eblg.getEBISL();
    const DisjointBoxLayout& dbl   = eblg.getDBL();
    const DataIterator&      dit   = dbl.dataIterator();
    const int                nbox  = dit.size();

    LayoutData<IntVectSet> irregCells(dbl);

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      irregCells[din] = ebisl[din].getIrregIVS(dbl[din]);
    }

    m_sigma[lvl] = RefCountedPtr<LevelData<BaseIVFAB<Real>>>(
      new LevelData<BaseIVFAB<Real>>(dbl, 1, IntVect::Zero, BaseIVFactory<Real>(ebisl, irregCells)));

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      (*m_sigma[lvl])[dit[mybox]].setVal(0.0);
    }
  }

  // The interface cells of a box: cut cells in both phases that are rows in both.
  auto getInterfaceCells = [&](const int a_level, const DataIndex& a_din) -> IntVectSet {
    IntVectSet interfaceCells;

    if (hasTwoPhases) {
      const Box                    box        = mflg[a_level]->getGrids()[a_din];
      const BaseFab<PetscAMRCell>& cells      = (*m_petscGrid->getAMRToPetsc()[a_level])[a_din];
      const IntVectSet             candidates = mflg[a_level]->interfaceRegion(box, a_din);

      for (IVSIterator ivsIt(candidates); ivsIt.ok(); ++ivsIt) {
        if (cells(ivsIt()).getPetscRow(0) >= 0 && cells(ivsIt()).getPetscRow(1) >= 0) {
          interfaceCells |= ivsIt();
        }
      }
    }

    return interfaceCells;
  };

  // The rows that the stencils of the local interface cells reach. They can be on other ranks, so their values are
  // gathered into a local vector.
  std::vector<PetscInt> rows;

  for (int lvl = 0; lvl < numLevels; lvl++) {
    const DisjointBoxLayout& dbl  = mflg[lvl]->getGrids();
    const DataIterator&      dit  = dbl.dataIterator();
    const int                nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const IntVectSet interfaceCells = getInterfaceCells(lvl, din);

      std::vector<PetscInt> boxRows;

      for (IVSIterator ivsIt(interfaceCells); ivsIt.ok(); ++ivsIt) {
        for (int k = 0; k < 2; k++) {
          Real         w;
          PetscStencil S;

          this->getNormalDerivativeStencil(w, S, k, lvl, din, VolIndex(ivsIt(), 0));

          for (int i = 0; i < S.size(); i++) {
            boxRows.push_back(S.row(i));
          }
        }
      }

#pragma omp critical(PetscHelmholtzSaturationChargeJumpBC_rows)
      {
        rows.insert(rows.end(), boxRows.begin(), boxRows.end());
      }
    }
  }

  std::sort(rows.begin(), rows.end());
  rows.erase(std::unique(rows.begin(), rows.end()), rows.end());

  // Collective: every rank takes part, also one without interface cells.
  IS         rowSet;
  Vec        localPhi;
  VecScatter scatter;

  PetscCallVoid(ISCreateGeneral(PETSC_COMM_SELF, rows.size(), rows.data(), PETSC_COPY_VALUES, &rowSet));
  PetscCallVoid(VecCreateSeq(PETSC_COMM_SELF, rows.size(), &localPhi));
  PetscCallVoid(VecScatterCreate(a_phi, rowSet, localPhi, nullptr, &scatter));
  PetscCallVoid(VecScatterBegin(scatter, a_phi, localPhi, INSERT_VALUES, SCATTER_FORWARD));
  PetscCallVoid(VecScatterEnd(scatter, a_phi, localPhi, INSERT_VALUES, SCATTER_FORWARD));

  const PetscScalar* phiValues = nullptr;

  PetscCallVoid(VecGetArrayRead(localPhi, &phiValues));

  // Apply a stencil to the gathered values.
  auto apply = [&](const PetscStencil& a_stencil) -> Real {
    Real ret = 0.0;

    for (int i = 0; i < a_stencil.size(); i++) {
      const auto idx = std::lower_bound(rows.begin(), rows.end(), a_stencil.row(i)) - rows.begin();

      ret += a_stencil.weight(i) * phiValues[idx];
    }

    return ret;
  };

  // sigma = b_0 D_0 + b_1 D_1 = b_o D_o, with phi_B = -S_z/w_z from D_z = 0.
  for (int lvl = 0; lvl < numLevels; lvl++) {
    const DisjointBoxLayout& dbl  = mflg[lvl]->getGrids();
    const DataIterator&      dit  = dbl.dataIterator();
    const int                nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const IntVectSet interfaceCells = getInterfaceCells(lvl, din);

      for (IVSIterator ivsIt(interfaceCells); ivsIt.ok(); ++ivsIt) {
        const VolIndex vof(ivsIt(), 0);

        Real         w[2];
        PetscStencil S[2];

        for (int k = 0; k < 2; k++) {
          this->getNormalDerivativeStencil(w[k], S[k], k, lvl, din, vof);
        }

        const Real bco = (*a_bCoefIrreg[lvl])[din].getIVFAB(o)(vof, 0);
        // Without a stencil on side z there is no flux through the interface, and no surface charge (see
        // getInterfaceFluxStencil).
        Real sigma = 0.0;

        if (w[z] > 0.0) {
          const Real phiB = -apply(S[z]) / w[z];

          sigma = bco * (w[o] * phiB + apply(S[o]));
        }

        (*m_sigma[lvl])[din](vof, 0) = sigma;
      }
    }
  }

  PetscCallVoid(VecRestoreArrayRead(localPhi, &phiValues));
  PetscCallVoid(VecScatterDestroy(&scatter));
  PetscCallVoid(VecDestroy(&localPhi));
  PetscCallVoid(ISDestroy(&rowSet));
}

const EBAMRIVData&
PetscHelmholtzSaturationChargeJumpBC::getSigma() const noexcept
{
  CH_TIME("PetscHelmholtzSaturationChargeJumpBC::getSigma");
  if (m_verbose) {
    pout() << "PetscHelmholtzSaturationChargeJumpBC::getSigma" << endl;
  }

  return m_sigma;
}

#include <CD_NamespaceFooter.H>

#endif
