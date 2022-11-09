/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBMGRestrict.cpp
  @brief  Implementation of CD_EBMGRestrict.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBCellFactory.H>
#include <CH_Timer.H>

// Our includes
#include <CD_BoxLoops.H>
#include <CD_DataOps.H>
#include <CD_EBMGRestrict.H>
#include <CD_NamespaceHeader.H>

EBMGRestrict::EBMGRestrict() noexcept
{
  CH_TIME("EBMGRestrict::EBMGRestrict(default)");

  m_isDefined = false;
}

EBMGRestrict::EBMGRestrict(const EBLevelGrid&   a_eblgFine,
                           const EBLevelGrid&   a_eblgCoar,
                           const ProblemDomain& a_domainCoar,
                           const int&           a_refRat) noexcept
{
  CH_TIME("EBMGRestrict::EBMGRestrict(full)");

  this->define(a_eblgFine, a_eblgCoar, a_domainCoar, a_refRat);
}

EBMGRestrict::~EBMGRestrict() noexcept { CH_TIME("EBMGRestrict::~EBMGRestrict"); }

void
EBMGRestrict::define(const EBLevelGrid&   a_eblgFine,
                     const EBLevelGrid&   a_eblgCoar,
                     const ProblemDomain& a_domainCoar,
                     const int&           a_refRat) noexcept
{
  CH_TIME("EBMGRestrict::define");

  m_refRat   = a_refRat;
  m_eblgFine = a_eblgFine;
  m_eblgCoar = a_eblgCoar;

#if 1 // Debug hook
  if (!a_eblgFine.getDBL().coarsenable(m_refRat)) {
    MayDay::Abort("EBMGRestrict::define -- the input grid is not coarsenable. We need to adopt a new strategy!");
  }
#endif

  // Create the coarsend layout
  coarsen(m_eblgCoFi, m_eblgFine, m_refRat);
  m_eblgCoFi.setMaxRefinementRatio(m_refRat);

  // Define buffer data.
  const DisjointBoxLayout& dblCoFi   = m_eblgCoFi.getDBL();
  const EBISLayout&        ebislCoFi = m_eblgCoFi.getEBISL();

  m_dataCoFi.define(dblCoFi, 1, IntVect::Zero, EBCellFactory(ebislCoFi));

  // Figure out the EB restriction stencils.
  m_vofitCoar.define(dblCoFi);
  m_restrictStencils.define(dblCoFi);

  for (DataIterator dit(dblCoFi); dit.ok(); ++dit) {
    const Box&       cellBox = dblCoFi[dit()];
    const EBISBox&   ebisBox = ebislCoFi[dit()];
    const EBGraph&   ebgraph = ebisBox.getEBGraph();
    const IntVectSet irreg   = ebisBox.getIrregIVS(cellBox);

    BaseIVFAB<VoFStencil>& restrictStencils = m_restrictStencils[dit()];
    VoFIterator&           vofit            = m_vofitCoar[dit()];

    vofit.define(irreg, ebgraph);
    restrictStencils.define(irreg, ebgraph, 1);

    for (vofit.reset(); vofit.ok(); ++vofit) {
      const VolIndex&        coarVoF     = vofit();
      const Vector<VolIndex> fineVoFs    = ebislCoFi.refine(coarVoF, m_refRat, dit());
      const int              numFineVoFs = 1. / fineVoFs.size();

      CH_assert(numFineVoFs > 0);

      VoFStencil& restrictSten = restrictStencils(coarVoF, 0);

      restrictSten.clear();
      for (int i = 0; i < numFineVoFs; i++) {
        restrictSten.add(fineVoFs[i], 1. / numFineVoFs);
      }
    }
  }
}

void
EBMGRestrict::restrict(LevelData<EBCellFAB>&       a_coarData,
                       const LevelData<EBCellFAB>& a_fineData,
                       const Interval              a_variables) const noexcept
{
  CH_TIME("EBMGRestrict::restrict");

  CH_assert(a_coarData.nComp() > a_variables.end());
  CH_assert(a_fineData.nComp() > a_variables.end());

  // Needed for single-valued data kernels.
  const Box refineBox(IntVect::Zero, (m_refRat - 1) * IntVect::Unit);
  const int weight = 1.0 / refineBox.numPts();

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {

    DataOps::setValue(m_dataCoFi, 0.0);

    const DisjointBoxLayout& dblCoar = m_eblgCoFi.getDBL();

    for (DataIterator dit(dblCoar); dit.ok(); ++dit) {
      EBCellFAB&       coarData = m_dataCoFi[dit()];
      const EBCellFAB& fineData = a_fineData[dit()];

      FArrayBox&       coarDataReg = coarData.getFArrayBox();
      const FArrayBox& fineDataReg = fineData.getFArrayBox();

      const BaseIVFAB<VoFStencil>& restrictStencils = m_restrictStencils[dit()];

      // Regular kernel.
      auto regularKernel = [&](const IntVect& ivCoar) -> void {
        for (BoxIterator bit(refineBox); bit.ok(); ++bit) {
          const IntVect ivFine = m_refRat * ivCoar + bit();

          coarDataReg(ivCoar, 0) += weight * fineDataReg(ivFine, ivar);
        }
      };

      // Irregular kernel
      auto irregularKernel = [&](const VolIndex& coarVoF) -> void {
        const VoFStencil& restrictSten = restrictStencils(coarVoF, 0);

        coarData(coarVoF, 0) = 0.0;
        for (int i = 0; i < restrictSten.size(); i++) {
          const VolIndex& fineVoF    = restrictSten.vof(i);
          const Real&     fineWeight = restrictSten.weight(i);

          coarData(coarVoF, 0) += fineWeight * fineData(fineVoF, ivar);
        }
      };

      // Run kernels
      const Box    coarBox  = dblCoar[dit()];
      VoFIterator& coarVoFs = m_vofitCoar[dit()];

      BoxLoops::loop(coarBox, regularKernel);
      BoxLoops::loop(coarVoFs, irregularKernel);
    }

    const Interval srcComps = Interval(0, 0);
    const Interval dstComps = Interval(ivar, ivar);

    m_dataCoFi.copyTo(srcComps, a_coarData, dstComps);
  }
}

#include <CD_NamespaceFooter.H>
