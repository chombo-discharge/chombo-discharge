/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBMGProlong.cpp
  @brief  Implementation of CD_EBMGProlong.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBCellFactory.H>
#include <CH_Timer.H>

// Our includes
#include <CD_BoxLoops.H>
#include <CD_DataOps.H>
#include <CD_EBMGProlong.H>
#include <CD_NamespaceHeader.H>

EBMGProlong::EBMGProlong() noexcept
{
  CH_TIME("EBMGProlong::EBMGProlong(default)");

  m_isDefined = false;
}

EBMGProlong::EBMGProlong(const EBLevelGrid& a_eblgFine, const EBLevelGrid& a_eblgCoar, const int& a_refRat) noexcept
{
  CH_TIME("EBMGProlong::EBMGProlong(full)");

  m_isDefined = false;

  this->define(a_eblgFine, a_eblgCoar, a_refRat);
}

EBMGProlong::~EBMGProlong() noexcept { CH_TIME("EBMGProlong::~EBMGProlong"); }

void
EBMGProlong::define(const EBLevelGrid& a_eblgFine, const EBLevelGrid& a_eblgCoar, const int& a_refRat) noexcept
{
  CH_TIME("EBMGProlong::define");

  CH_assert(a_refRat >= 2);
  CH_assert(a_refRat % 2 == 0);

  m_refRat   = a_refRat;
  m_eblgFine = a_eblgFine;
  m_eblgCoar = a_eblgCoar;

  if (!(a_eblgFine.getDBL().coarsenable(m_refRat))) {
    MayDay::Abort("EBMGProlong::define -- the input grid is not coarsenable. We need to adopt a new strategy!");
  }

  // Create the coarsened layout
  coarsen(m_eblgCoFi, m_eblgFine, m_refRat);
  m_eblgCoFi.setMaxRefinementRatio(m_refRat);

  const DisjointBoxLayout& dblCoFi = m_eblgCoFi.getDBL();
  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();

  const EBISLayout& ebislCoFi = m_eblgCoFi.getEBISL();
  const EBISLayout& ebislFine = m_eblgFine.getEBISL();

  // Figure out the EB restriction stencils.
  m_vofitFine.define(dblCoFi);
  m_prolongStencils.define(dblCoFi);

  for (DataIterator dit(dblFine); dit.ok(); ++dit) {
    const Box&       cellBox = dblFine[dit()];
    const EBISBox&   ebisBox = ebislFine[dit()];
    const EBGraph&   ebgraph = ebisBox.getEBGraph();
    const IntVectSet irreg   = ebisBox.getIrregIVS(cellBox);

    BaseIVFAB<VoFStencil>& prolongStencils = m_prolongStencils[dit()];
    VoFIterator&           vofit           = m_vofitFine[dit()];

    vofit.define(irreg, ebgraph);
    prolongStencils.define(irreg, ebgraph, 1);

    for (vofit.reset(); vofit.ok(); ++vofit) {
      const VolIndex& fineVoF = vofit();
      const VolIndex  coarVoF = ebislFine.coarsen(fineVoF, m_refRat, dit());

      VoFStencil& prolongSten = prolongStencils(fineVoF, 0);

      prolongSten.clear();
      prolongSten.add(coarVoF, 1.0);
    }
  }

  // Note: MUST have the same number of ghost cells as the buffer being defined in prolongResidual
  m_copier.ghostDefine(m_eblgCoar.getDBL(), dblCoFi, m_eblgCoar.getDomain(), IntVect::Zero);

  m_isDefined = true;
}

void
EBMGProlong::prolongResidual(LevelData<EBCellFAB>&       a_fineData,
                             const LevelData<EBCellFAB>& a_coarData,
                             const Interval              a_variables) const noexcept
{
  CH_TIMERS("EBMGProlong::prolongResidual");
  CH_TIMER("EBMGProlong::prolongResidual::regular_cells", t1);
  CH_TIMER("EBMGProlong::prolongResidual::irregular_cells", t2);

  CH_assert(m_isDefined);
  CH_assert(a_fineData.nComp() > a_variables.end());
  CH_assert(a_coarData.nComp() > a_variables.end());

  LevelData<EBCellFAB> coFiData(m_eblgCoFi.getDBL(), 1, IntVect::Zero, EBCellFactory(m_eblgCoFi.getEBISL()));

  const Box refineBox(IntVect::Zero, (m_refRat - 1) * IntVect::Unit);

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {

    // Copy coarse data to buffer.
    const Interval srcComps = Interval(ivar, ivar);
    const Interval dstComps = Interval(0, 0);

    a_coarData.copyTo(srcComps, coFiData, dstComps, m_copier);

    // Add coarse-grid residual to the fine grid. Recall that m_eblgCoFi is a coarsening of the fine grid
    // so this runs over the part of the coarse level that is covered by the finer level.
    const DisjointBoxLayout& dblCoFi   = m_eblgCoFi.getDBL();
    const EBISLayout&        ebislFine = m_eblgFine.getEBISL();

    for (DataIterator dit(dblCoFi); dit.ok(); ++dit) {
      EBCellFAB&       fineData = a_fineData[dit()];
      const EBCellFAB& coarData = coFiData[dit()];

      FArrayBox&       fineDataReg = fineData.getFArrayBox();
      const FArrayBox& coarDataReg = coarData.getFArrayBox();

      const EBISBox& ebisBoxFine = ebislFine[dit()];

      const BaseIVFAB<VoFStencil>& prolongStencils = m_prolongStencils[dit()];

      // Regular kernel.
      auto regularKernel = [&](const IntVect& ivCoar) -> void {
        for (BoxIterator bit(refineBox); bit.ok(); ++bit) {
          const IntVect ivFine = m_refRat * ivCoar + bit();

          // Put a guard for cut-cells on the fine grid because the irregular will do those.
          if (!(ebisBoxFine.isIrregular(ivFine))) {
            fineDataReg(ivFine, ivar) += coarDataReg(ivCoar, 0);
          }
        }
      };

      // Irregular kernel
      auto irregularKernel = [&](const VolIndex& fineVoF) -> void {
        const VoFStencil& prolongSten = prolongStencils(fineVoF, 0);

        for (int i = 0; i < prolongSten.size(); i++) {
          const VolIndex& coarVoF    = prolongSten.vof(i);
          const Real&     coarWeight = prolongSten.weight(i);

          fineData(fineVoF, ivar) += coarWeight * coarData(coarVoF, 0);
        }
      };

      // Run kernels
      const Box    coarBox  = dblCoFi[dit()];
      VoFIterator& fineVoFs = m_vofitFine[dit()];

      CH_START(t1);
      BoxLoops::loop(coarBox, regularKernel);
      CH_STOP(t1);

      CH_START(t2);
      BoxLoops::loop(fineVoFs, irregularKernel);
      CH_STOP(t2);
    }
  }
}

#include <CD_NamespaceFooter.H>
