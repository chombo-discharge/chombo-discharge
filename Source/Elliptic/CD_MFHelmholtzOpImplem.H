/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzOpImplem.H
  @brief  Implementation of CD_MFHelmholtzOp.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFHelmholtzOp.H>
#include <CD_ParallelOps.H>
#include <CD_NamespaceHeader.H>

template <typename Duration>
Vector<long long>
MFHelmholtzOp::computeOperatorLoads(LevelData<MFCellFAB>& a_phi, const int a_numApply) noexcept
{
  CH_TIME("MFHelmholtzOp::computeOperatorLoads");

  // TLDR: This routine estimates the time spent in each grid patch for a typical relaxation step. This includes
  //       coarse-fine interpolation, BC matching, and applying the operator.

  CH_assert(a_numApply > 0);
  CH_assert(a_phi.ghostVect() == m_ghostPhi);

  LevelData<MFCellFAB> Lphi;
  this->create(Lphi, a_phi);

  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();
  const DataIterator&      dit = dbl.dataIterator();

  Vector<long long> loads(dbl.size(), 0LL);

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din     = dit[mybox];
    const int        intCode = din.intCode();

    const auto t0 = std::chrono::steady_clock::now();

    for (int i = 0; i < a_numApply; i++) {

      // Homogeneous interpolation with coarser level.
      if (m_hasCoar) {
        for (auto& op : m_helmOps) {
          const int iphase = op.first;

          RefCountedPtr<EBMultigridInterpolator>& phaseInterpolator = m_interpolator.getInterpolator(iphase);

          EBCellFAB& phi = (EBCellFAB&)a_phi[din].getPhase(iphase);

          phaseInterpolator->coarseFineInterpH(phi, Interval(m_comp, m_comp), din);
        }
      }

      // Matching time
      if (m_multifluid) {
        m_jumpBC->matchBC((*m_jump)[din], a_phi[din], true, din);
      }

      // Apply operator application
      for (auto& op : m_helmOps) {
        const Box cellBox = Lphi.disjointBoxLayout()[din];

        const int iphase = op.first;

        EBCellFAB& Lph = Lphi[din].getPhase(iphase);
        EBCellFAB& phi = a_phi[din].getPhase(iphase);

        const EBCellFAB&       Acoef      = (*m_Acoef)[din].getPhase(iphase);
        const EBFluxFAB&       Bcoef      = (*m_Bcoef)[din].getPhase(iphase);
        const BaseIVFAB<Real>& BcoefIrreg = *(*m_BcoefIrreg)[din].getPhasePtr(iphase);

        op.second->applyOp(Lph, phi, Acoef, Bcoef, BcoefIrreg, cellBox, din, true);
      }
    }

    const auto t1       = std::chrono::steady_clock::now();
    const auto duration = (std::chrono::duration_cast<Duration>(t1 - t0));

    loads[intCode] = (long long)duration.count();
  }

  ParallelOps::vectorSum(loads);

  return loads;
}

#include <CD_NamespaceFooter.H>
