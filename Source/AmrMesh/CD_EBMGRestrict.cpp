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

EBMGRestrict::EBMGRestrict(const EBLevelGrid& a_eblgFine, const EBLevelGrid& a_eblgCoar, const int& a_refRat) noexcept
{
  CH_TIME("EBMGRestrict::EBMGRestrict(full)");

  m_isDefined = false;

  this->define(a_eblgFine, a_eblgCoar, a_refRat);
}

EBMGRestrict::~EBMGRestrict() noexcept
{
  CH_TIME("EBMGRestrict::~EBMGRestrict");
}

void
EBMGRestrict::define(const EBLevelGrid& a_eblgFine, const EBLevelGrid& a_eblgCoar, const int& a_refRat) noexcept
{
  CH_TIME("EBMGRestrict::define");

  CH_assert(a_refRat >= 2);
  CH_assert(a_refRat % 2 == 0);

  m_refRat   = a_refRat;
  m_eblgFine = a_eblgFine;
  m_eblgCoar = a_eblgCoar;

  if (!(a_eblgFine.getDBL().coarsenable(m_refRat))) {
    MayDay::Abort("EBMGRestrict::define -- the input grid is not coarsenable. We need to adopt a new strategy!");
  }

  // Create the coarsend layout
  coarsen(m_eblgCoFi, m_eblgFine, m_refRat);
  m_eblgCoFi.setMaxRefinementRatio(m_refRat);

  // Define buffer data.
  const DisjointBoxLayout& dblCoFi   = m_eblgCoFi.getDBL();
  const EBISLayout&        ebislCoFi = m_eblgCoFi.getEBISL();

  // Weighting factor for stencils.
  const Real fineWeight = 1. / std::pow(m_refRat, SpaceDim);

  // Figure out the EB restriction stencils.
  m_vofitCoar.define(dblCoFi);
  m_restrictStencils.define(dblCoFi);

  // Note: MUST have the same number of ghost cells as the buffer being defined in the restriction function.
  m_copier.define(dblCoFi, m_eblgCoar.getDBL(), IntVect::Zero);

  const DataIterator& dit  = dblCoFi.dataIterator();
  const int           nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const Box&       cellBox = dblCoFi[din];
    const EBISBox&   ebisBox = ebislCoFi[din];
    const EBGraph&   ebgraph = ebisBox.getEBGraph();
    const IntVectSet irreg   = ebisBox.getIrregIVS(cellBox);

    BaseIVFAB<VoFStencil>& restrictStencils = m_restrictStencils[din];
    VoFIterator&           vofit            = m_vofitCoar[din];

    vofit.define(irreg, ebgraph);
    restrictStencils.define(irreg, ebgraph, 1);

    for (vofit.reset(); vofit.ok(); ++vofit) {
      const VolIndex&        coarVoF     = vofit();
      const Vector<VolIndex> fineVoFs    = ebislCoFi.refine(coarVoF, m_refRat, din);
      const int              numFineVoFs = fineVoFs.size();

      VoFStencil& restrictSten = restrictStencils(coarVoF, 0);

      restrictSten.clear();
      for (int i = 0; i < numFineVoFs; i++) {
        restrictSten.add(fineVoFs[i], fineWeight);
      }
    }
  }

  m_isDefined = true;
}

void
EBMGRestrict::restrictResidual(LevelData<EBCellFAB>&       a_coarData,
                               const LevelData<EBCellFAB>& a_fineData,
                               const Interval              a_variables) const noexcept
{
  CH_TIMERS("EBMGRestrict::restrictResidual");
  CH_TIMER("EBMGRestrict::restrictResidual::regular_cells", t2);
  CH_TIMER("EBMGRestrict::restrictResidual::irregular_cells", t3);

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_variables.end());
  CH_assert(a_fineData.nComp() > a_variables.end());

  // Needed for single-valued data kernels.
  const Box  refineBox     = Box(IntVect::Zero, (m_refRat - 1) * IntVect::Unit);
  const Real regFineWeight = 1.0 / std::pow(m_refRat, SpaceDim);

  const DisjointBoxLayout& dblCoFi   = m_eblgCoFi.getDBL();
  const EBISLayout&        ebislCoFi = m_eblgCoFi.getEBISL();

  LevelData<EBCellFAB> coFiData(dblCoFi, 1, IntVect::Zero, EBCellFactory(ebislCoFi));

  const DataIterator& dit  = dblCoFi.dataIterator();
  const int           nbox = dit.size();

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din      = dit[mybox];
      EBCellFAB&       coarData = coFiData[din];
      const EBCellFAB& fineData = a_fineData[din];

      FArrayBox&       coarDataReg = coarData.getFArrayBox();
      const FArrayBox& fineDataReg = fineData.getFArrayBox();

      const BaseIVFAB<VoFStencil>& restrictStencils = m_restrictStencils[din];

      // Regular kernel.
      auto regularKernel = [&](const IntVect& ivCoar) -> void {
        for (BoxIterator bit(refineBox); bit.ok(); ++bit) {
          const IntVect ivFine = m_refRat * ivCoar + bit();

          coarDataReg(ivCoar, 0) += regFineWeight * fineDataReg(ivFine, ivar);
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
      coarData.setVal(0.0);

      const Box    coarBox  = dblCoFi[din];
      VoFIterator& coarVoFs = m_vofitCoar[din];

      CH_START(t2);
      BoxLoops::loop(coarBox, regularKernel);
      CH_STOP(t2);

      CH_START(t3);
      BoxLoops::loop(coarVoFs, irregularKernel);
      CH_STOP(t3);
    }

    const Interval srcComps = Interval(0, 0);
    const Interval dstComps = Interval(ivar, ivar);

    coFiData.copyTo(srcComps, a_coarData, dstComps, m_copier);
  }
}

#include <CD_NamespaceFooter.H>
