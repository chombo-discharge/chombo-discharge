/* chombo-discharge
 * Copyright © 2025 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBNonConservativeDivergence.cpp
  @brief  Implementation of CD_EBNonConservativeDiverge.cpp
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBNonConservativeDivergence.H>
#include <CD_BoxLoops.H>
#include <CD_VofUtils.H>
#include <CD_NamespaceHeader.H>

EBNonConservativeDivergence::EBNonConservativeDivergence() noexcept
{
  CH_TIME("EBNonConservativeDivergence::EBNonConservativeDivergence(weak)");

  m_isDefined = false;
}

EBNonConservativeDivergence::EBNonConservativeDivergence(const EBLevelGrid& a_eblg, const int a_radius) noexcept
{
  CH_TIME("EBNonConservativeDivergence::EBNonConservativeDivergence(strong)");

  CH_assert(a_eblg.isDefined());
  CH_assert(a_radius >= 1);

  m_isDefined = false;

  this->define(a_eblg, a_radius);
}

EBNonConservativeDivergence::~EBNonConservativeDivergence() noexcept
{
  CH_TIME("EBNonConservativeDivergence::~EBNonConservativeDivergence");
}

void
EBNonConservativeDivergence::define(const EBLevelGrid& a_eblg, const int a_radius) noexcept
{
  CH_TIME("EBNonConservativeDivergence::define");

  CH_assert(a_eblg.isDefined());
  CH_assert(a_radius >= 1);

  m_eblg = a_eblg;

  const DisjointBoxLayout& dbl    = m_eblg.getDBL();
  const EBISLayout&        ebisl  = m_eblg.getEBISL();
  const ProblemDomain&     domain = m_eblg.getDomain();
  const DataIterator&      dit    = dbl.dataIterator();

  m_vofIterator.define(dbl);
  m_stencils.define(dbl);

  const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const EBISBox&   ebisBox  = ebisl[din];
    const Box&       cellBox  = dbl[din];
    const EBGraph&   ebGraph  = ebisBox.getEBGraph();
    const IntVectSet irregIVS = ebisBox.getIrregIVS(cellBox);

    VoFIterator&           vofit    = m_vofIterator[din];
    BaseIVFAB<VoFStencil>& stencils = m_stencils[din];

    vofit.define(irregIVS, ebGraph);
    stencils.define(irregIVS, ebGraph, 1);

    for (vofit.reset(); vofit.ok(); ++vofit) {
      const VolIndex& vof     = vofit();
      VoFStencil&     stencil = stencils(vof, 0);

      Real sumKappa = 0.0;

      const Vector<VolIndex> neighborVofs = VofUtils::getVofsInRadius(vof,
                                                                      ebisBox,
                                                                      a_radius,
                                                                      VofUtils::Connectivity::MonotonePath,
                                                                      true);

      for (int i = 0; i < neighborVofs.size(); i++) {
        const VolIndex& ivof  = neighborVofs[i];
        const Real&     kappa = ebisBox.volFrac(ivof);

        sumKappa += kappa;

        stencil.add(ivof, kappa);
      }

      stencil *= 1. / sumKappa;
    }
  }

  m_isDefined = true;
}

void
EBNonConservativeDivergence::nonConservativeDivergence(LevelData<BaseIVFAB<Real>>& a_nonConsDivF,
                                                       const LevelData<EBCellFAB>& a_kappaDivF) const noexcept
{
  CH_TIME("EBNonConservativeDivergence::nonConservativeDivergence");

  CH_assert(m_isDefined);
  CH_assert(a_nonConsDivF.isDefined());
  CH_assert(a_kappaDivF.isDefined());
  CH_assert(a_nonConsDivF.nComp() == a_kappaDivF.nComp());
  CH_assert(a_nonConsDivF.nComp() >= 1);
  CH_assert(a_nonConsDivF.disjointBoxLayout() == m_eblg.getDBL());
  CH_assert(a_kappaDivF.disjointBoxLayout() == m_eblg.getDBL());

  const DisjointBoxLayout& dbl = m_eblg.getDBL();
  const DataIterator&      dit = dbl.dataIterator();

  const int nComp = a_nonConsDivF.nComp();
  const int nbox  = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex&             din      = dit[mybox];
    const BaseIVFAB<VoFStencil>& stencils = m_stencils[din];

    BaseIVFAB<Real>& nonConsDivF = a_nonConsDivF[din];
    const EBCellFAB& kappaDivF   = a_kappaDivF[din];

    for (int comp = 0; comp < nComp; comp++) {

      auto kernel = [&](const VolIndex& vof) -> void {
        nonConsDivF(vof, comp) = 0.0;

        const VoFStencil& stencil = stencils(vof, 0);
        for (int i = 0; i < stencil.size(); i++) {
          const VolIndex& ivof    = stencil.vof(i);
          const Real&     iweight = stencil.weight(i);

          nonConsDivF(vof, comp) += iweight * kappaDivF(ivof, comp);
        }
      };

      BoxLoops::loop(m_vofIterator[din], kernel);
    }
  }
}

#include <CD_NamespaceFooter.H>
