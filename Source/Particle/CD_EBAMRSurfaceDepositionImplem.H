/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_EBAMRSurfaceDepositionImplem.H
  @brief  Implementation of CD_EBAMRSurfaceDeposition.H
  @author Robert Marskar
*/

#ifndef CD_EBAMRSURFACEDEPOSITIONIMPLEM_H
#define CD_EBAMRSURFACEDEPOSITIONIMPLEM_H

// Std includes
#include <type_traits>

// Chombo includes
#include <ParmParse.H>
#include <CH_Timer.H>

// Our includes
#include <CD_IrregAddOp.H>
#include <CD_DataOps.H>
#include <CD_ParticleOps.H>
#include <CD_EBAMRSurfaceDeposition.H>
#include <CD_NamespaceHeader.H>

template <typename P, typename Traits>
void
EBAMRSurfaceDeposition::deposit(EBAMRIVData&                           a_meshData,
                                const ParticleContainerSoA<P, Traits>& a_particles) const noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::deposit(SoA)");
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::deposit(ParticleContainerSoA)" << endl;
  }

  CH_assert(a_meshData.getRealm() == a_particles.getRealm());

  // Deposit on this level
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = m_ebGrids[lvl]->getDBL();
    const EBISLayout&        ebisl = m_ebGrids[lvl]->getEBISL();
    const DataIterator&      dit   = dbl.dataIterator();

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box                    cellBox  = dbl[din];
      const EBISBox&               ebisbox  = ebisl[din];
      const BaseIVFAB<VoFStencil>& stencils = (*m_depositionStencils[lvl])[din];

      BaseIVFAB<Real>&              meshData = (*m_data[lvl])[din];
      const ParticleSoA<P, Traits>& leaf     = a_particles[lvl][din];

      meshData.setVal(0.0);

      for (std::size_t i = 0; i < leaf.size(); i++) {
        const IntVect iv = ParticleOps::getParticleCellIndex(leaf.position(i), m_probLo, m_dx[lvl]);

        // If this fails the particle is not in a valid cell.
        if (!(cellBox.contains(iv))) {
          MayDay::Error("CD_EBAMRSurfaceDeposition::deposit(SoA) -- particle is not inside the box");
        }

        if (ebisbox.isIrregular(iv)) {
          const VoFStencil& stencil = stencils(VolIndex(iv, 0), 0);

          for (int k = 0; k < stencil.size(); k++) {
            const VolIndex& stencilVoF       = stencil.vof(k);
            const Real&     stencilWeight    = stencil.weight(k);
            const Real      depositionWeight = stencilWeight * leaf.weight(i);

            meshData(stencilVoF, 0) += depositionWeight;
          }
        }
      }
    }

    // Above, we will have deposited over patch boundaries. Add the deposted data into the neighboring patches.
    m_data[lvl]->exchange(Interval(0, 0), m_copierLevel[lvl], IrregAddOp());
  }

  // Ensure conservation across coarse-fine interface. This involves interpolation of the data going from
  // the coarse level to the fine level, and coarsening of the data from the fine level to coarse level.
  this->addInvalidCoarseDataToFineData();
  this->addFineGhostDataToValidCoarData();

  // Finally, copy our buffers to the input data.
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    CH_assert(!(a_meshData[lvl].isNull()));
    CH_assert(a_meshData[lvl]->nComp() == 1);

    const Interval interv = Interval(0, 0);

    m_data[lvl]->copyTo(interv, *a_meshData[lvl], interv, m_validToValidCopiers[lvl]);
  }
}

#include <CD_NamespaceFooter.H>

#endif
