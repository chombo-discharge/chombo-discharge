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
#include <BaseIVFactory.H>

// Our includes
#include <CD_EBFineInterp.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

EBFineInterp::EBFineInterp()
{
  CH_TIME("EBFineInterp::EBFineInterp(weak)");
  m_isDefined = false;
}

EBFineInterp::EBFineInterp(const EBLevelGrid&        a_eblgFine,
                           const EBLevelGrid&        a_eblgCoar,
                           const int&                a_refRat,
                           const int&                a_nComp,
                           const EBIndexSpace* const a_ebisPtr)
{
  CH_TIME("EBFineInterp::EBFineInterp(full)");

  this->define(a_eblgFine, a_eblgCoar, a_refRat, a_nComp, a_ebisPtr);
}

EBFineInterp::~EBFineInterp() {}

void
EBFineInterp::define(const EBLevelGrid&        a_eblgFine,
                     const EBLevelGrid&        a_eblgCoar,
                     const int&                a_refRat,
                     const int&                a_nComp,
                     const EBIndexSpace* const a_ebisPtr)
{
  CH_TIME("EBFineInterp::define");

  CH_assert(a_refRat > 1);
  CH_assert(a_nComp > 0);
  CH_assert(!(a_ebisPtr == nullptr));

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

  // Define the irreg data holder.
  LayoutData<IntVectSet> irregCells;

  for (DataIterator dit(m_coarsenedFineGrids); dit.ok(); ++dit) {
    const Box&     cellBox = m_coarsenedFineGrids[dit()];
    const EBISBox& ebisbox = m_coarsenedFineEBISL[dit()];

    irregCells[dit()] = ebisbox.getIrregIVS(cellBox);
  }

  // Should not need ghost cells for this one (because the irreg regrids don't use slopes but keeps it local).
  m_coarsenedFineIrregData.define(m_coarsenedFineGrids,
                                  1,
                                  IntVect::Zero,
                                  BaseIVFactory<Real>(m_coarsenedFineEBISL, irregCells));
}

void
EBFineInterp::regridNoSlopes(LevelData<EBCellFAB>&       a_fineData,
                             const LevelData<EBCellFAB>& a_coarData,
                             const Interval&             a_variables)
{
  CH_TIME("EBFineInterp::regridNoSlopes");

  CH_assert(m_isDefined);

  a_coarData.copyTo(a_variables, m_coarsenedFineData, a_variables);

  for (DataIterator dit(m_coarsenedFineGrids); dit.ok(); ++dit) {
    this->regridNoSlopes(a_fineData[dit()], m_coarsenedFineData[dit()], dit(), a_variables);
  }
}

void
EBFineInterp::regridNoSlopes(EBCellFAB&       a_fineData,
                             const EBCellFAB& a_coarData,
                             const DataIndex& a_dit,
                             const Interval&  a_variables)
{
  CH_TIME("EBFineInterp::regridNoSlopes");

  const Box coarBox = m_coarsenedFineGrids[a_dit];
  const Box refiBox = Box(IntVect::Zero, (m_refRat - 1) * IntVect::Unit);

  FArrayBox&       fineDataReg = a_fineData.getFArrayBox();
  const FArrayBox& coarDataReg = a_coarData.getFArrayBox();

  const EBISBox& ebisBoxCoar = m_coarsenedFineEBISL[a_dit];
  const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[a_dit];

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {

    // Regular kernel. Set the fine data equal to the coarse data.
    auto regularKernel = [&](const IntVect& coarIV) -> void {
      const Real& coarVal = coarDataReg(coarIV, ivar);

      for (BoxIterator bit(refiBox); bit.ok(); ++bit) {
        const IntVect& fineIV     = m_refRat * coarIV + bit();
        fineDataReg(fineIV, ivar) = coarVal;
      }
    };

    auto irregularKernel = [&](const VolIndex& coarVoF) -> void {
      const Vector<VolIndex> fineVoFs = m_coarsenedFineEBISL.refine(coarVoF, m_refRat, a_dit);

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

      // Initialize irregular data.
      for (int i = 0; i < fineVoFs.size(); i++) {
        a_fineData(fineVoFs[i], ivar) = kappaFactor * a_coarData(coarVoF, ivar);
      }
    };

    // Irregular kernel region.
    VoFIterator vofit(m_irregRegions[a_dit], ebisBoxCoar.getEBGraph());

    BoxLoops::loop(coarBox, regularKernel);
    BoxLoops::loop(vofit, irregularKernel);
  }
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

  MayDay::Abort("EBFineInterp::regridConservative -- not implemented");
}

void
EBFineInterp::regridArithmetic(LevelData<BaseIVFAB<Real>>&       a_fineData,
                               const LevelData<BaseIVFAB<Real>>& a_coarData,
                               const Interval&                   a_variables)
{
  CH_TIME("EBFineInterp::regridArithmetic");

  MayDay::Abort("EBFineInterp::regridArithmetic -- not implemented");
}

#include <CD_NamespaceFooter.H>
