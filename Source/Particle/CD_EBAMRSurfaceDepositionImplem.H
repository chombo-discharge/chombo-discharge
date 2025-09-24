/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBAMRSurfaceDepositionImplem.H
  @brief  Implementation of CD_EBAMRSurfaceDeposition.H
  @author Robert Marskar
*/

#ifndef CD_EBAMRSurfaceDepositionImplem_H
#define CD_EBAMRSurfaceDepositionImplem_H

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

template <class P, class Ret, Ret (P::*MemberFunc)() const>
void
EBAMRSurfaceDeposition::deposit(EBAMRIVData& a_meshData, const ParticleContainer<P>& a_particles) const noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::deposit");
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::deposit<P, const Real& P::*func const>" << endl;
  }

  using NoRef = typename std::remove_reference<Ret>::type;
  using Base  = typename std::remove_cv<NoRef>::type;

  static_assert(std::is_same<Base, Real>::value, "Ret should be Real or const Real&");

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

      BaseIVFAB<Real>& meshData  = (*m_data[lvl])[din];
      const List<P>&   particles = a_particles[lvl][din].listItems();

      meshData.setVal(0.0);

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        const P& p = lit();

        const IntVect iv = ParticleOps::getParticleCellIndex(p.position(), m_probLo, m_dx[lvl]);

        // If this fails the particle is not in a valid cell.
        if (!(cellBox.contains(iv))) {
          MayDay::Error("CD_EBAMRSurfaceDeposition::deposit -- particle is not inside the box");
        }

        if (ebisbox.isIrregular(iv)) {
          const VoFStencil& stencil = stencils(VolIndex(iv, 0), 0);

          for (int i = 0; i < stencil.size(); i++) {
            const VolIndex& stencilVoF       = stencil.vof(i);
            const Real&     stencilWeight    = stencil.weight(i);
            const Real      depositionWeight = stencilWeight * (p.*MemberFunc)();

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
