/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBFineInterp.cpp
  @brief  Implementation of EBFineInterp.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <EBCellFactory.H>
#include <BaseIVFactory.H>

// Our includes
#include <CD_EBFineInterp.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

EBFineInterp::EBFineInterp() noexcept
{
  CH_TIME("EBFineInterp::EBFineInterp(weak)");
  
  m_isDefined = false;
}

EBFineInterp::EBFineInterp(const EBLevelGrid&        a_eblgFine,
                           const EBLevelGrid&        a_eblgCoFi,
                           const EBLevelGrid&        a_eblgCoar,
                           const int&                a_refRat,
                           const int&                a_nComp,
                           const EBIndexSpace* const a_ebisPtr) noexcept
{
  CH_TIME("EBFineInterp::EBFineInterp(full)");

  this->define(a_eblgFine, a_eblgCoFi, a_eblgCoar, a_refRat, a_nComp, a_ebisPtr);
}

EBFineInterp::~EBFineInterp() noexcept {}

void
EBFineInterp::define(const EBLevelGrid&        a_eblgFine,
                     const EBLevelGrid&        a_eblgCoFi,		     
                     const EBLevelGrid&        a_eblgCoar,
                     const int&                a_refRat,
                     const int&                a_nComp,
                     const EBIndexSpace* const a_ebisPtr) noexcept
{
  CH_TIME("EBFineInterp::define");

  CH_assert(a_refRat > 1);
  CH_assert(a_nComp > 0);
  CH_assert(a_ebisPtr != nullptr);

  EBPWLFineInterp::define(a_eblgFine.getDBL(),
                          a_eblgCoar.getDBL(),
                          a_eblgFine.getEBISL(),
                          a_eblgCoar.getEBISL(),
                          a_eblgCoar.getDomain(),
                          a_refRat,
                          a_nComp,
                          a_ebisPtr);

  m_eblgFine = a_eblgFine;
  m_eblgCoar = a_eblgCoar;
  m_eblgCoFi = a_eblgCoFi;

  // Define the irreg data holder.
  LayoutData<IntVectSet> irregFine(m_eblgFine.getDBL());
  LayoutData<IntVectSet> irregCoar(m_coarsenedFineGrids);

  m_fineVoFs.define(m_eblgFine.getDBL());

  for (DataIterator dit(m_eblgFine.getDBL()); dit.ok(); ++dit) {
    const Box& fineBox = m_eblgFine.getDBL()[dit()];
    const Box& coarBox = m_coarsenedFineGrids[dit()];

    const EBISBox& fineEBISBox = m_eblgFine.getEBISL()[dit()];
    const EBISBox& coarEBISBox = m_coarsenedFineEBISL[dit()];

    irregFine[dit()] = fineEBISBox.getIrregIVS(fineBox);
    irregCoar[dit()] = coarEBISBox.getIrregIVS(coarBox);

    m_fineVoFs[dit()].define(irregFine[dit()], fineEBISBox.getEBGraph());
  }

  // Should not need ghost cells for this one (because the irreg regrids don't use slopes but keeps it local).
  m_irregCoFi.define(m_coarsenedFineGrids, 1, IntVect::Zero, BaseIVFactory<Real>(m_coarsenedFineEBISL, irregCoar));

  m_conservativeWeights.define(m_eblgFine.getDBL(),
                               1,
                               IntVect::Zero,
                               BaseIVFactory<Real>(m_eblgFine.getEBISL(), irregFine));

  this->defineWeights();
}

void
EBFineInterp::defineWeights() noexcept
{
  CH_TIME("EBFineInterp::defineWeights");

  const EBISLayout& ebislFine = m_eblgFine.getEBISL();
  const EBISLayout& ebislCoar = m_coarsenedFineEBISL;

  const Real areaFactor = std::pow(m_refRat, SpaceDim - 1);

  for (DataIterator dit(m_eblgFine.getDBL()); dit.ok(); ++dit) {
    const EBISBox& ebisBoxFine = ebislFine[dit()];
    const EBISBox& ebisBoxCoar = ebislCoar[dit()];

    BaseIVFAB<Real>& weights = m_conservativeWeights[dit()];

    auto kernel = [&](const VolIndex& fineVoF) -> void {
      // TLDR: We want refRat^(D-1) * coarArea/fineArea. We refine the coarse vof (which also gives regular cells)
      // but that won't matter because they don't have a boundary area.

      const VolIndex         coarVoF     = ebislFine.coarsen(fineVoF, m_refRat, dit());
      const Vector<VolIndex> refCoarVoFs = ebislCoar.refine(coarVoF, m_refRat, dit());

      Real fineArea = 0.0;
      Real coarArea = ebisBoxCoar.bndryArea(coarVoF);

      for (int i = 0; i < refCoarVoFs.size(); i++) {
        if (ebisBoxFine.isIrregular(refCoarVoFs[i].gridIndex())) {
          fineArea += ebisBoxFine.bndryArea(refCoarVoFs[i]);
        }
      }

      if (fineArea > 0.0) {
        weights(fineVoF, 0) = areaFactor * coarArea / fineArea;
      }
      else {
        weights(fineVoF, 0) = 0.0;
      }
    };

    BoxLoops::loop(m_fineVoFs[dit()], kernel);
  }
}

void
EBFineInterp::regridConservativeNoSlopes(LevelData<EBCellFAB>&       a_fineData,
                                         const LevelData<EBCellFAB>& a_coarData,
                                         const Interval&             a_variables) const noexcept
{
  CH_TIMERS("EBFineInterp::regridConservativeNoSlopes(LD<EBCellFAB>)");
  CH_TIMER("EBFineInterp::define_buffer", t1);
  CH_TIMER("EBFineInterp::copyTo", t2);

  pout() << "in here" << endl;

  CH_assert(m_isDefined);
  CH_assert(a_fineData.nComp() > a_variables.end());
  CH_assert(a_coarData.nComp() > a_variables.end());

  CH_START(t1);
  const DisjointBoxLayout& dblCoFi   = m_eblgCoFi.getDBL();
  const EBISLayout&        ebislCoFi = m_eblgCoFi.getEBISL();
  LevelData<EBCellFAB>     coarsenedFineData(dblCoFi, 1, IntVect::Zero, EBCellFactory(ebislCoFi));
  CH_STOP(t1);

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {
    CH_START(t2);
    a_coarData.copyTo(Interval(ivar, ivar), coarsenedFineData, Interval(0, 0));
    CH_STOP(t2);

    for (DataIterator dit(dblCoFi); dit.ok(); ++dit) {
      this->regridConservativeNoSlopes(a_fineData[dit()], coarsenedFineData[dit()], dit(), ivar, 0);
    }
  }
}

void
EBFineInterp::regridConservativeNoSlopes(EBCellFAB&       a_fineData,
                                         const EBCellFAB& a_coarData,
                                         const DataIndex& a_dit,
                                         const int&       a_fineVar,
                                         const int&       a_coarVar) const noexcept
{
  CH_TIMERS("EBFineInterp::regridConservativeNoSlopes(EBCellFAB)");
  CH_TIMER("EBFineInterp::regular_regrid", t1);
  CH_TIMER("EBFineInterp::irregular_regrid", t2);

  CH_assert(a_fineData.nComp() > a_fineVar);
  CH_assert(a_coarData.nComp() > a_coarVar);

  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();
  const DisjointBoxLayout& dblCoar = m_eblgCoFi.getDBL();

  const EBISLayout& ebislFine = m_eblgFine.getEBISL();
  const EBISLayout& ebislCoar = m_eblgCoFi.getEBISL();

  const EBISBox& ebisBoxFine = ebislFine[a_dit];
  const EBISBox& ebisBoxCoar = ebislCoar[a_dit];

  const Box coarBox = dblCoar[a_dit];
  const Box refiBox = Box(IntVect::Zero, (m_refRat - 1) * IntVect::Unit);

  const Real volFactor = std::pow(m_refRat, SpaceDim);

  FArrayBox&       fineDataReg = a_fineData.getFArrayBox();
  const FArrayBox& coarDataReg = a_coarData.getFArrayBox();

  // Regular kernel. Set the fine data equal to the coarse data.
  auto regularKernel = [&](const IntVect& coarIV) -> void {
    const Real& coarVal = coarDataReg(coarIV, a_coarVar);

    for (BoxIterator bit(refiBox); bit.ok(); ++bit) {
      const IntVect& fineIV = m_refRat * coarIV + bit();

      fineDataReg(fineIV, a_fineVar) = coarVal;
    }
  };

  // Cut-cell kernel.
  auto irregularKernel = [&](const VolIndex& coarVoF) -> void {
    const Vector<VolIndex> fineVoFs = ebislCoar.refine(coarVoF, m_refRat, a_dit);

    // Compute kappaCoar/sum(kappaFine)
    Real kappaFactor = 0.0;
    for (int i = 0; i < fineVoFs.size(); i++) {
      kappaFactor += ebisBoxFine.volFrac(fineVoFs[i]);
    }

    if (kappaFactor > 1.E-9) {
      kappaFactor = ebisBoxCoar.volFrac(coarVoF) / kappaFactor;
    }
    else {
      kappaFactor = 1.0;
    }

    kappaFactor *= volFactor;

    // Initialize irregular data.
    for (int i = 0; i < fineVoFs.size(); i++) {
      a_fineData(fineVoFs[i], a_fineVar) = kappaFactor * a_coarData(coarVoF, a_coarVar);
    }
  };

  CH_START(t1);
  BoxLoops::loop(coarBox, regularKernel);
  CH_STOP(t1);

  CH_START(t2);
  VoFIterator vofit(m_irregRegions[a_dit], ebisBoxCoar.getEBGraph());
  BoxLoops::loop(vofit, irregularKernel);
  CH_STOP(t2);
}

void
EBFineInterp::regridMinMod(LevelData<EBCellFAB>&       a_fineData,
                           const LevelData<EBCellFAB>& a_coarData,
                           const Interval&             a_variables)
{
  CH_TIME("EBFineInterp::regridMinMod");

  EBPWLFineInterp::interpolate(a_fineData, a_coarData, a_variables);
}

void
EBFineInterp::regridConservative(LevelData<BaseIVFAB<Real>>&       a_fineData,
                                 const LevelData<BaseIVFAB<Real>>& a_coarData,
                                 const Interval&                   a_variables)
{
  CH_TIME("EBFineInterp::regridConservative");

  CH_assert(a_fineData.nComp() >= a_variables.size());
  CH_assert(a_coarData.nComp() >= a_variables.size());

  CH_assert(a_fineData.nComp() > a_variables.end());
  CH_assert(a_coarData.nComp() > a_variables.end());

  CH_assert(a_fineData.disjointBoxLayout() == m_eblgFine.getDBL());
  CH_assert(a_coarData.disjointBoxLayout() == m_eblgCoar.getDBL());

  for (int comp = a_variables.begin(); comp <= a_variables.end(); comp++) {
    a_coarData.copyTo(Interval(comp, comp), m_irregCoFi, Interval(0, 0));

    const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();

    const EBISLayout& ebislFine = m_eblgFine.getEBISL();
    const EBISLayout& ebislCoar = m_coarsenedFineEBISL;

    for (DataIterator dit(dblFine); dit.ok(); ++dit) {
      const EBISBox& fineEBISBox = ebislFine[dit()];
      const EBISBox& coarEBISBox = ebislCoar[dit()];

      BaseIVFAB<Real>&       fineData = a_fineData[dit()];
      const BaseIVFAB<Real>& coarData = m_irregCoFi[dit()];
      const BaseIVFAB<Real>& weights  = m_conservativeWeights[dit()];

      auto kernel = [&](const VolIndex& fineVoF) -> void {
        const VolIndex coarVoF = ebislFine.coarsen(fineVoF, m_refRat, dit());

        CH_assert(coarEBISBox.isIrregular(coarVoF.gridIndex()));
        CH_assert(fineEBISBox.isIrregular(fineVoF.gridIndex()));

        fineData(fineVoF, comp) = coarData(coarVoF, 0) * weights(fineVoF, 0);
      };

      BoxLoops::loop(m_fineVoFs[dit()], kernel);
    }
  }
}

void
EBFineInterp::regridArithmetic(LevelData<BaseIVFAB<Real>>&       a_fineData,
                               const LevelData<BaseIVFAB<Real>>& a_coarData,
                               const Interval&                   a_variables)
{
  CH_TIME("EBFineInterp::regridArithmetic");

  CH_assert(a_fineData.nComp() >= a_variables.size());
  CH_assert(a_coarData.nComp() >= a_variables.size());

  CH_assert(a_fineData.nComp() > a_variables.end());
  CH_assert(a_coarData.nComp() > a_variables.end());

  CH_assert(a_fineData.disjointBoxLayout() == m_eblgFine.getDBL());
  CH_assert(a_coarData.disjointBoxLayout() == m_eblgCoar.getDBL());

  for (int comp = a_variables.begin(); comp <= a_variables.end(); comp++) {
    a_coarData.copyTo(Interval(comp, comp), m_irregCoFi, Interval(0, 0));

    const EBISLayout& ebislFine = m_eblgFine.getEBISL();

    for (DataIterator dit(m_eblgFine.getDBL()); dit.ok(); ++dit) {
      BaseIVFAB<Real>&       fineData = a_fineData[dit()];
      const BaseIVFAB<Real>& coarData = m_irregCoFi[dit()];

      auto kernel = [&](const VolIndex& fineVoF) -> void {
        const VolIndex coarVoF = ebislFine.coarsen(fineVoF, m_refRat, dit());

        fineData(fineVoF, comp) = coarData(coarVoF, 0);
      };

      BoxLoops::loop(m_fineVoFs[dit()], kernel);
    }
  }
}

#include <CD_NamespaceFooter.H>
