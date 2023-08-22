/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBCoarAve.cpp
  @brief  Implementation of CD_EBCoarAve.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBCellFactory.H>
#include <EBFluxFactory.H>
#include <BaseIVFactory.H>
#include <CH_Timer.H>

// Our includes
#include <CD_EBCoarAve.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

EBCoarAve::EBCoarAve() noexcept
{
  CH_TIME("EBCoarAve::EBCoarAve");
  m_isDefined = false;
}

EBCoarAve::~EBCoarAve() noexcept
{
  CH_TIME("EBCoarAve::~EBCoarAve");

  m_cellCopiers.clear();
  m_faceCopiers.clear();
  m_ebCopiers.clear();
}

EBCoarAve::EBCoarAve(const DisjointBoxLayout& a_dblFine,
                     const DisjointBoxLayout& a_dblCoar,
                     const EBISLayout&        a_ebislFine,
                     const EBISLayout&        a_ebislCoar,
                     const ProblemDomain&     a_domainCoar,
                     const int&               a_refRat,
                     const EBIndexSpace*      a_ebisPtr) noexcept
{
  CH_TIME("EBCoarAve::EBCoarAve(DBL version)");

  CH_assert(a_ebisPtr->isDefined());
  CH_assert(a_refRat >= 2);
  CH_assert(a_refRat % 2 == 0);

  ProblemDomain domainFine = a_domainCoar;
  domainFine.refine(a_refRat);

  EBLevelGrid eblgFine(a_dblFine, a_ebislFine, domainFine);
  EBLevelGrid eblgCoar(a_dblCoar, a_ebislCoar, a_domainCoar);
  EBLevelGrid eblgCoFi;

  coarsen(eblgCoFi, eblgFine, a_refRat);
  eblgCoFi.getEBISL().setMaxRefinementRatio(a_refRat, a_ebisPtr);

  this->define(eblgFine, eblgCoar, eblgCoFi, a_refRat);
}

EBCoarAve::EBCoarAve(const EBLevelGrid& a_eblgFine,
                     const EBLevelGrid& a_eblgCoar,
                     const EBLevelGrid& a_eblgCoFi,
                     const int&         a_refRat) noexcept
{
  CH_TIME("EBCoarAve::EBCoarAve(EBLevelGrid version)");

  m_isDefined = false;

  this->define(a_eblgFine, a_eblgCoar, a_eblgCoFi, a_refRat);
}

void
EBCoarAve::define(const EBLevelGrid& a_eblgFine,
                  const EBLevelGrid& a_eblgCoar,
                  const EBLevelGrid& a_eblgCoFi,
                  const int&         a_refRat) noexcept
{
  CH_TIME("EBCoarAve::define");

  CH_assert(a_refRat >= 2);
  CH_assert(a_refRat % 2 == 0);

  m_cellCopiers.clear();
  m_faceCopiers.clear();
  m_ebCopiers.clear();

  m_eblgFine = a_eblgFine;
  m_eblgCoar = a_eblgCoar;
  m_eblgCoFi = a_eblgCoFi;
  m_refRat   = a_refRat;

  m_irregSetsCoFi.define(m_eblgCoFi.getDBL());

  const DisjointBoxLayout& dblCoFi = m_eblgCoFi.getDBL();
  const DataIterator&      ditCoFi = dblCoFi.dataIterator();

  const int nbox = ditCoFi.size();
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = ditCoFi[mybox];

    m_irregSetsCoFi[din] = m_eblgCoFi.getEBISL()[din].getIrregIVS(dblCoFi[din]);
  }

  this->defineCellStencils();
  this->defineFaceStencils();
  this->defineEBStencils();
  this->defineBuffers();

  m_isDefined = true;
}

void
EBCoarAve::defineCellStencils() noexcept
{
  CH_TIME("EBCoarAve::defineCellStencils");

  const DisjointBoxLayout& dblCoar = m_eblgCoFi.getDBL();
  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();

  const DataIterator& ditCoar = dblCoar.dataIterator();

  const EBISLayout& ebislCoar = m_eblgCoFi.getEBISL();
  const EBISLayout& ebislFine = m_eblgFine.getEBISL();

  const Real dxCoar   = 1.0;
  const Real dxFine   = dxCoar / m_refRat;
  const Real dxFactor = std::pow(dxFine / dxCoar, SpaceDim);
  const int  nbox     = ditCoar.size();

  m_irregCellsCoFi.define(dblCoar);
  m_cellConservativeStencils.define(dblCoar);
  m_cellArithmeticStencils.define(dblCoar);
  m_cellHarmonicStencils.define(dblCoar);

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din         = ditCoar[mybox];
    const EBISBox&   ebisBoxCoar = ebislCoar[din];
    const EBISBox&   ebisBoxFine = ebislFine[din];
    const EBGraph&   ebGraphCoar = ebisBoxCoar.getEBGraph();

    VoFIterator&           irregCells           = m_irregCellsCoFi[din];
    BaseIVFAB<VoFStencil>& conservativeStencils = m_cellConservativeStencils[din];
    BaseIVFAB<VoFStencil>& arithmeticStencils   = m_cellArithmeticStencils[din];
    BaseIVFAB<VoFStencil>& harmonicStencils     = m_cellHarmonicStencils[din];

    irregCells.define(m_irregSetsCoFi[din], ebGraphCoar);
    conservativeStencils.define(m_irregSetsCoFi[din], ebGraphCoar, 1);
    arithmeticStencils.define(m_irregSetsCoFi[din], ebGraphCoar, 1);
    harmonicStencils.define(m_irregSetsCoFi[din], ebGraphCoar, 1);

    auto buildStencils = [&](const VolIndex& coarVoF) -> void {
      const Real             kappaC      = ebisBoxCoar.volFrac(coarVoF);
      const Vector<VolIndex> fineVoFs    = ebislCoar.refine(coarVoF, m_refRat, din);
      const int              numFineVoFs = fineVoFs.size();

      VoFStencil& arithSten = arithmeticStencils(coarVoF, 0);
      VoFStencil& harmSten  = harmonicStencils(coarVoF, 0);
      VoFStencil& consSten  = conservativeStencils(coarVoF, 0);

      arithSten.clear();
      harmSten.clear();
      consSten.clear();

      if (numFineVoFs > 0) {

        for (int ifine = 0; ifine < numFineVoFs; ifine++) {
          const VolIndex& fineVoF = fineVoFs[ifine];
          const Real      kappaF  = ebisBoxFine.volFrac(fineVoF);

          arithSten.add(fineVoF, 1.0 / numFineVoFs);
          harmSten.add(fineVoF, 1.0 / numFineVoFs);

          if (kappaC > 0.0) {
            consSten.add(fineVoF, kappaF * dxFactor / kappaC);
          }
          else {
            consSten.add(fineVoF, 1.0 / numFineVoFs);
          }
        }
      }
    };

    BoxLoops::loop(irregCells, buildStencils);
  }
}

void
EBCoarAve::defineFaceStencils() noexcept
{
  CH_TIME("EBCoarAve::defineFaceStencils");

  const DisjointBoxLayout& dblCoar = m_eblgCoFi.getDBL();
  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();

  const DataIterator& ditCoar = dblCoar.dataIterator();

  const EBISLayout& ebislCoar = m_eblgCoFi.getEBISL();
  const EBISLayout& ebislFine = m_eblgFine.getEBISL();

  const Real dxCoar   = 1.0;
  const Real dxFine   = dxCoar / m_refRat;
  const Real dxFactor = std::pow(dxFine / dxCoar, SpaceDim - 1);
  const int  nbox     = ditCoar.size();

  m_irregFacesCoFi.define(dblCoar);
  m_faceArithmeticStencils.define(dblCoar);
  m_faceHarmonicStencils.define(dblCoar);
  m_faceConservativeStencils.define(dblCoar);

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = ditCoar[mybox];

    const EBISBox& ebisBoxCoar = ebislCoar[din];
    const EBISBox& ebisBoxFine = ebislFine[din];
    const EBGraph& ebGraphCoar = ebisBoxCoar.getEBGraph();

    const IntVectSet irregIVS = ebisBoxCoar.getIrregIVS(dblCoar[din]);

    for (int dir = 0; dir < SpaceDim; dir++) {
      FaceIterator& faceIt = m_irregFacesCoFi[din][dir];

      BaseIFFAB<FaceStencil>& arithmeticStencils   = m_faceArithmeticStencils[din][dir];
      BaseIFFAB<FaceStencil>& harmonicStencils     = m_faceHarmonicStencils[din][dir];
      BaseIFFAB<FaceStencil>& conservativeStencils = m_faceConservativeStencils[din][dir];

      faceIt.define(irregIVS, ebGraphCoar, dir, FaceStop::SurroundingWithBoundary);

      conservativeStencils.define(irregIVS, ebGraphCoar, dir, 1);
      arithmeticStencils.define(irregIVS, ebGraphCoar, dir, 1);
      harmonicStencils.define(irregIVS, ebGraphCoar, dir, 1);

      auto buildStencils = [&](const FaceIndex& coarFace) -> void {
        const Real              areaCoar     = ebisBoxCoar.areaFrac(coarFace);
        const Vector<FaceIndex> fineFaces    = ebislCoar.refine(coarFace, m_refRat, din);
        const int               numFineFaces = fineFaces.size();

        FaceStencil& arithSten = arithmeticStencils(coarFace, 0);
        FaceStencil& harmSten  = harmonicStencils(coarFace, 0);
        FaceStencil& consSten  = conservativeStencils(coarFace, 0);

        arithSten.clear();
        harmSten.clear();
        consSten.clear();

        if (numFineFaces > 0) {
          for (int ifine = 0; ifine < numFineFaces; ifine++) {
            const FaceIndex& fineFace = fineFaces[ifine];
            const Real       areaFine = ebisBoxFine.areaFrac(fineFace);

            arithSten.add(fineFace, 1.0 / numFineFaces);
            harmSten.add(fineFace, 1.0 / numFineFaces);

            if (areaCoar > 0.0) {
              consSten.add(fineFace, areaFine * dxFactor / areaCoar);
            }
            else {
              consSten.add(fineFace, 1.0 / numFineFaces);
            }
          }
        }
      };

      BoxLoops::loop(faceIt, buildStencils);
    }
  }
}

void
EBCoarAve::defineEBStencils() noexcept
{
  CH_TIME("EBCoarAve::defineEBStencils");

  const DisjointBoxLayout& dblCoar = m_eblgCoFi.getDBL();
  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();

  const DataIterator& ditCoar = dblCoar.dataIterator();

  const EBISLayout& ebislCoar = m_eblgCoFi.getEBISL();
  const EBISLayout& ebislFine = m_eblgFine.getEBISL();

  const Real dxCoar   = 1.0;
  const Real dxFine   = dxCoar / m_refRat;
  const Real dxFactor = std::pow(dxFine / dxCoar, SpaceDim - 1);
  const int  nbox     = ditCoar.size();

  m_ebArithmeticStencils.define(dblCoar);
  m_ebHarmonicStencils.define(dblCoar);
  m_ebConservativeStencils.define(dblCoar);

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = ditCoar[mybox];

    const EBISBox& ebisBoxCoar = ebislCoar[din];
    const EBISBox& ebisBoxFine = ebislFine[din];
    const EBGraph& ebGraphCoar = ebisBoxCoar.getEBGraph();

    BaseIVFAB<VoFStencil>& conservativeStencils = m_ebConservativeStencils[din];
    BaseIVFAB<VoFStencil>& arithmeticStencils   = m_ebArithmeticStencils[din];
    BaseIVFAB<VoFStencil>& harmonicStencils     = m_ebHarmonicStencils[din];

    conservativeStencils.define(m_irregSetsCoFi[din], ebGraphCoar, 1);
    arithmeticStencils.define(m_irregSetsCoFi[din], ebGraphCoar, 1);
    harmonicStencils.define(m_irregSetsCoFi[din], ebGraphCoar, 1);

    auto buildStencils = [&](const VolIndex& coarVoF) -> void {
      const Real             areaCoar    = ebisBoxCoar.bndryArea(coarVoF);
      const Vector<VolIndex> refinedVoFs = ebislCoar.refine(coarVoF, m_refRat, din);

      VoFStencil& arithSten = arithmeticStencils(coarVoF, 0);
      VoFStencil& harmSten  = harmonicStencils(coarVoF, 0);
      VoFStencil& consSten  = conservativeStencils(coarVoF, 0);

      arithSten.clear();
      harmSten.clear();
      consSten.clear();

      // Get the irregular cells coming from coarsening of coarVoF
      Vector<VolIndex> fineVoFs;
      for (int i = 0; i < refinedVoFs.size(); i++) {
        if (ebisBoxFine.isIrregular(refinedVoFs[i].gridIndex())) {
          fineVoFs.push_back(refinedVoFs[i]);
        }
      }
      const int numFineVoFs = fineVoFs.size();

      if (numFineVoFs > 0) {
        for (int ifine = 0; ifine < numFineVoFs; ifine++) {
          const VolIndex& fineVoF  = fineVoFs[ifine];
          const Real      areaFine = ebisBoxFine.bndryArea(fineVoF);

          arithSten.add(fineVoF, 1.0 / numFineVoFs);
          harmSten.add(fineVoF, 1.0 / numFineVoFs);

          if (areaCoar > 0.0) {
            consSten.add(fineVoF, areaFine * dxFactor / areaCoar);
          }
          else {
            consSten.add(fineVoF, 1.0 / numFineVoFs);
          }
        }
      }
    };

    VoFIterator& irregCells = m_irregCellsCoFi[din];

    BoxLoops::loop(irregCells, buildStencils);
  }
}

void
EBCoarAve::defineBuffers() noexcept
{
  CH_TIME("EBCoarAve::defineBuffers");
}

void
EBCoarAve::averageData(LevelData<EBCellFAB>&       a_coarData,
                       const LevelData<EBCellFAB>& a_fineData,
                       const Interval&             a_variables,
                       const Average&              a_average) const noexcept
{
  CH_TIMERS("EBCoarAve::averageData(LD<EBCellFAB>)");
  CH_TIMER("EBCoarAve::averageData(LD<EBCellFAB>)::average", t1);
  CH_TIMER("EBCoarAve::averageData(LD<EBCellFAB>)::copyTo", t2);

  CH_assert(m_isDefined);
  CH_assert(a_fineData.nComp() > a_variables.end());
  CH_assert(a_coarData.nComp() > a_variables.end());

  const DisjointBoxLayout& dblCoFi = m_eblgCoFi.getDBL();
  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();

  const DataIterator& dit  = dblFine.dataIterator();
  const int           nbox = dit.size();

  LevelData<EBCellFAB> coFiData(dblCoFi, 1, IntVect::Zero, EBCellFactory(m_eblgCoFi.getEBISL()));

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {

    CH_START(t1);
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB&       coarData = coFiData[din];
      const EBCellFAB& fineData = a_fineData[din];

      // Switch between methods.
      switch (a_average) {
      case Average::Arithmetic: {
        this->arithmeticAverage(coarData, fineData, din, 0, ivar);

        break;
      }
      case Average::Harmonic: {
        this->harmonicAverage(coarData, fineData, din, 0, ivar);

        break;
      }
      case Average::Conservative: {
        this->conservativeAverage(coarData, fineData, din, 0, ivar);

        break;
      }
      default: {
        MayDay::Error("EBCoarAve::averageData(LD<EBCellFAB>) - logic bust");

        break;
      }
      }
    }
    CH_STOP(t1);

    CH_START(t2);
    const Interval srcInterv = Interval(0, 0);
    const Interval dstInterv = Interval(ivar, ivar);

    Copier& cellCopier = m_cellCopiers[a_coarData.ghostVect()];
    if (!cellCopier.isDefined()) {
      cellCopier.define(m_eblgCoFi.getDBL(), m_eblgCoar.getDBL(), a_coarData.ghostVect());
    }

    coFiData.copyTo(srcInterv, a_coarData, dstInterv, cellCopier);
    CH_STOP(t2);
  }
}

void
EBCoarAve::arithmeticAverage(EBCellFAB&       a_coarData,
                             const EBCellFAB& a_fineData,
                             const DataIndex& a_datInd,
                             const int&       a_coarVar,
                             const int&       a_fineVar) const noexcept
{
  CH_TIMERS("EBCoarAve::arithmeticAverage(EBCellFAB)");
  CH_TIMER("EBCoarAve::arithmeticAverage(EBCellFAB)::regular_cells", t1);
  CH_TIMER("EBCoarAve::arithmeticAverage(EBCellFAB)::irregular_cells", t2);

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_coarVar);
  CH_assert(a_fineData.nComp() > a_fineVar);

  const Real dxCoar   = 1.0;
  const Real dxFine   = dxCoar / m_refRat;
  const Real dxFactor = std::pow(dxFine / dxCoar, SpaceDim);
  const Box  refiBox  = Box(IntVect::Zero, (m_refRat - 1) * IntVect::Unit);

  FArrayBox&       coarDataReg = a_coarData.getFArrayBox();
  const FArrayBox& fineDataReg = a_fineData.getFArrayBox();

  const BaseIVFAB<VoFStencil>& coarseningStencils = m_cellArithmeticStencils[a_datInd];

  // Coarsening of regular cells.
  auto regularKernel = [&](const IntVect& iv) -> void {
    coarDataReg(iv, a_coarVar) = 0.0;

    for (BoxIterator bit(refiBox); bit.ok(); ++bit) {
      coarDataReg(iv, a_coarVar) += fineDataReg(m_refRat * iv + bit(), a_fineVar);
    }

    coarDataReg(iv, a_coarVar) *= dxFactor;
  };

  // Coarsening of cut-cells.
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    a_coarData(vof, a_coarVar) = 0.0;

    const VoFStencil& stencil = coarseningStencils(vof, 0);

    for (int i = 0; i < stencil.size(); i++) {
      a_coarData(vof, a_coarVar) += stencil.weight(i) * a_fineData(stencil.vof(i), a_fineVar);
    }
  };

  CH_START(t1);
  const Box& coarBox = m_eblgCoFi.getDBL()[a_datInd];

  BoxLoops::loop(coarBox, regularKernel);
  CH_STOP(t1);

  CH_START(t2);
  VoFIterator& vofit = m_irregCellsCoFi[a_datInd];

  BoxLoops::loop(vofit, irregularKernel);
  CH_STOP(t2);
}

void
EBCoarAve::harmonicAverage(EBCellFAB&       a_coarData,
                           const EBCellFAB& a_fineData,
                           const DataIndex& a_datInd,
                           const int&       a_coarVar,
                           const int&       a_fineVar) const noexcept

{
  CH_TIMERS("EBCoarAve::harmonicAverage(EBCellFAB)");
  CH_TIMER("EBCoarAve::harmonicAverage(EBCellFAB)::regular_cells", t1);
  CH_TIMER("EBCoarAve::harmonicAverage(EBCellFAB)::irregular_cells", t2);

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_coarVar);
  CH_assert(a_fineData.nComp() > a_fineVar);

  // Regular cells
  const Box  refiBox    = Box(IntVect::Zero, (m_refRat - 1) * IntVect::Unit);
  const Real numPerCoar = refiBox.numPts();

  FArrayBox&       coarDataReg = a_coarData.getFArrayBox();
  const FArrayBox& fineDataReg = a_fineData.getFArrayBox();

  const BaseIVFAB<VoFStencil>& coarseningStencils = m_cellHarmonicStencils[a_datInd];

  // Harmonic coarsening of regular cells. Harmonic averaging is phiCoar = n/sum_{i<n}(1/x_i)
  auto regularKernel = [&](const IntVect& iv) -> void {
    Real coarVal = 0.0;

    for (BoxIterator bit(refiBox); bit.ok(); ++bit) {
      coarVal += 1.0 / fineDataReg(m_refRat * iv + bit(), a_fineVar);
    }

    coarDataReg(iv, a_coarVar) = numPerCoar / coarVal;
  };

  // Coarsening of cut-cells. We've put the stencil weights = 1/n so we only need to accumulate and invert.
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    a_coarData(vof, a_coarVar) = 0.0;

    const VoFStencil& stencil = coarseningStencils(vof, 0);

    for (int i = 0; i < stencil.size(); i++) {
      a_coarData(vof, a_coarVar) += stencil.weight(i) / a_fineData(stencil.vof(i), a_fineVar);
    }

    a_coarData(vof, a_coarVar) = 1. / a_coarData(vof, a_coarVar);
  };

  CH_START(t1);
  const Box& coarBox = m_eblgCoFi.getDBL()[a_datInd];

  BoxLoops::loop(coarBox, regularKernel);
  CH_STOP(t1);

  // Irregular cells
  CH_START(t2);
  VoFIterator& vofit = m_irregCellsCoFi[a_datInd];

  BoxLoops::loop(vofit, irregularKernel);
  CH_STOP(t2);
}

void
EBCoarAve::conservativeAverage(EBCellFAB&       a_coarData,
                               const EBCellFAB& a_fineData,
                               const DataIndex& a_datInd,
                               const int&       a_coarVar,
                               const int&       a_fineVar) const noexcept
{
  CH_TIMERS("EBCoarAve::conservativeAverage(EBCellFAB)");
  CH_TIMER("EBCoarAve::conservativeAverage(EBCellFAB)::regular_cells", t1);
  CH_TIMER("EBCoarAve::conservativeAverage(EBCellFAB)::irregular_cells", t2);

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_coarVar);
  CH_assert(a_fineData.nComp() > a_fineVar);

  const Real dxCoar   = 1.0;
  const Real dxFine   = dxCoar / m_refRat;
  const Real dxFactor = std::pow(dxFine / dxCoar, SpaceDim);
  const Box  refiBox  = Box(IntVect::Zero, (m_refRat - 1) * IntVect::Unit);

  FArrayBox&       coarDataReg = a_coarData.getFArrayBox();
  const FArrayBox& fineDataReg = a_fineData.getFArrayBox();

  const BaseIVFAB<VoFStencil>& coarseningStencils = m_cellConservativeStencils[a_datInd];

  // Coarsening of regular cells.
  auto regularKernel = [&](const IntVect& iv) -> void {
    coarDataReg(iv, a_coarVar) = 0.0;

    for (BoxIterator bit(refiBox); bit.ok(); ++bit) {
      coarDataReg(iv, a_coarVar) += fineDataReg(m_refRat * iv + bit(), a_fineVar);
    }

    coarDataReg(iv, a_coarVar) *= dxFactor;
  };

  // Coarsening of cut-cells.
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    a_coarData(vof, a_coarVar) = 0.0;

    const VoFStencil& stencil = coarseningStencils(vof, 0);

    for (int i = 0; i < stencil.size(); i++) {
      a_coarData(vof, a_coarVar) += stencil.weight(i) * a_fineData(stencil.vof(i), a_fineVar);
    }
  };

  CH_START(t1);
  const Box& coarBox = m_eblgCoFi.getDBL()[a_datInd];

  BoxLoops::loop(coarBox, regularKernel);
  CH_STOP(t1);

  CH_START(t2);
  VoFIterator& vofit = m_irregCellsCoFi[a_datInd];

  BoxLoops::loop(vofit, irregularKernel);
  CH_STOP(t2);
}

void
EBCoarAve::averageData(LevelData<EBFluxFAB>&       a_coarData,
                       const LevelData<EBFluxFAB>& a_fineData,
                       const Interval&             a_variables,
                       const Average&              a_average) const noexcept
{
  CH_TIMERS("EBCoarAve::averageData(LD<EBFluxFAB>)");
  CH_TIMER("EBCoarAve::averageData(LD<EBFluxFAB>)::average", t1);

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_variables.end());
  CH_assert(a_fineData.nComp() > a_variables.end());

  const DisjointBoxLayout& dblCoFi = m_eblgCoFi.getDBL();
  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();

  const DataIterator& dit  = dblFine.dataIterator();
  const int           nbox = dit.size();

  LevelData<EBFluxFAB> coFiData(dblCoFi, 1, IntVect::Zero, EBFluxFactory(m_eblgCoFi.getEBISL()));

  CH_START(t1);
  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      EBFluxFAB&       coarData = coFiData[din];
      const EBFluxFAB& fineData = a_fineData[din];

      switch (a_average) {
      case Average::Arithmetic: {
        for (int dir = 0; dir < SpaceDim; dir++) {
          this->arithmeticAverage(coarData[dir], fineData[dir], din, 0, ivar, dir);
        }

        break;
      }
      case Average::Harmonic: {
        for (int dir = 0; dir < SpaceDim; dir++) {
          this->harmonicAverage(coarData[dir], fineData[dir], din, 0, ivar, dir);
        }

        break;
      }
      case Average::Conservative: {
        for (int dir = 0; dir < SpaceDim; dir++) {
          this->conservativeAverage(coarData[dir], fineData[dir], din, 0, ivar, dir);
        }

        break;
      }
      default: {
        MayDay::Error("EBCoarAve::averageData(LD<EBFluxFAB>) - logic bust");

        break;
      }
      }
    }
    CH_STOP(t1);

    const Interval srcInterv = Interval(0, 0);
    const Interval dstInterv = Interval(ivar, ivar);

    Copier& faceCopier = m_faceCopiers[a_coarData.ghostVect()];
    if (!faceCopier.isDefined()) {
      faceCopier.define(m_eblgCoFi.getDBL(), m_eblgCoar.getDBL(), a_coarData.ghostVect());
    }

    coFiData.copyTo(srcInterv, a_coarData, dstInterv, faceCopier);
  }
}

void
EBCoarAve::arithmeticAverage(EBFaceFAB&       a_coarData,
                             const EBFaceFAB& a_fineData,
                             const DataIndex& a_datInd,
                             const int&       a_coarVar,
                             const int&       a_fineVar,
                             const int&       a_dir) const noexcept
{
  CH_TIMERS("EBCoarAve::arithmeticAverage(EBFaceFAB)");
  CH_TIMER("EBCoarAve::arithmeticAverage(EBFaceFAB)::regular_cells", t1);
  CH_TIMER("EBCoarAve::arithmeticAverage(EBFaceFAB)::irregular_cells", t2);

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_coarVar);
  CH_assert(a_fineData.nComp() > a_fineVar);

  FArrayBox&       coarDataReg = a_coarData.getFArrayBox();
  const FArrayBox& fineDataReg = a_fineData.getFArrayBox();

  const BaseIFFAB<FaceStencil>& coarseningStencils = m_faceArithmeticStencils[a_datInd].at(a_dir);

  // Kernel for doing regular faces.
  const Real invFinePerCoar = 1.0 / std::pow(m_refRat, SpaceDim - 1);
  const int  xDoLoop        = (a_dir == 0) ? 0 : 1;
  const int  yDoLoop        = (a_dir == 1) ? 0 : 1;
#if CH_SPACEDIM == 3
  const int zDoLoop = (a_dir == 2) ? 0 : 1;
#endif

  auto regularKernel = [&](const IntVect& iv) -> void {
    coarDataReg(iv, a_coarVar) = 0.0;

#if CH_SPACEDIM == 3
    for (int k = 0; k <= (m_refRat - 1) * zDoLoop; k++) {
#endif
      for (int j = 0; j <= (m_refRat - 1) * yDoLoop; j++) {
        for (int i = 0; i <= (m_refRat - 1) * xDoLoop; i++) {
          const IntVect ivFine = iv * m_refRat + IntVect(D_DECL(i, j, k));

          coarDataReg(iv, a_coarVar) += fineDataReg(ivFine, a_fineVar);
        }
      }
#if CH_SPACEDIM == 3
    }
#endif

    coarDataReg(iv, a_coarVar) *= invFinePerCoar;
  };

  // Kernel for doing irregular faces.
  auto irregularKernel = [&](const FaceIndex& face) -> void {
    a_coarData(face, a_coarVar) = 0.0;

    const FaceStencil& stencil = coarseningStencils(face, 0);
    for (int i = 0; i < stencil.size(); i++) {
      a_coarData(face, a_coarVar) += stencil.weight(i) * a_fineData(stencil.face(i), a_fineVar);
    }
  };

  // Regular faces.
  CH_START(t1);
  const Box& coarBox     = m_eblgCoFi.getDBL()[a_datInd];
  const Box  coarFaceBox = surroundingNodes(coarBox, a_dir);

  BoxLoops::loop(coarFaceBox, regularKernel);
  CH_STOP(t1);

  // Irregular faces
  CH_START(t2);
  FaceIterator& faceIt = m_irregFacesCoFi[a_datInd][a_dir];

  BoxLoops::loop(faceIt, irregularKernel);
  CH_STOP(t2);
}

void
EBCoarAve::harmonicAverage(EBFaceFAB&       a_coarData,
                           const EBFaceFAB& a_fineData,
                           const DataIndex& a_datInd,
                           const int&       a_coarVar,
                           const int&       a_fineVar,
                           const int&       a_dir) const noexcept
{
  CH_TIMERS("EBCoarAve::harmonicAverage(EBFaceFAB)");
  CH_TIMER("EBCoarAve::harmonicAverage(EBFaceFAB)::regular_cells", t1);
  CH_TIMER("EBCoarAve::harmonicAverage(EBFaceFAB)::irregular_cells", t2);

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_coarVar);
  CH_assert(a_fineData.nComp() > a_fineVar);

  FArrayBox&       coarDataReg = a_coarData.getFArrayBox();
  const FArrayBox& fineDataReg = a_fineData.getFArrayBox();

  const BaseIFFAB<FaceStencil>& coarseningStencils = m_faceHarmonicStencils[a_datInd].at(a_dir);

  const Real finePerCoar = std::pow(m_refRat, SpaceDim - 1);
  const int  xDoLoop     = (a_dir == 0) ? 0 : 1;
  const int  yDoLoop     = (a_dir == 1) ? 0 : 1;
#if CH_SPACEDIM == 3
  const int zDoLoop = (a_dir == 2) ? 0 : 1;
#endif

  // Kernel for doing regular faces.
  auto regularKernel = [&](const IntVect& iv) -> void {
    coarDataReg(iv, a_coarVar) = 0.0;

#if CH_SPACEDIM == 3
    for (int k = 0; k <= (m_refRat - 1) * zDoLoop; k++) {
#endif
      for (int j = 0; j <= (m_refRat - 1) * yDoLoop; j++) {
        for (int i = 0; i <= (m_refRat - 1) * xDoLoop; i++) {
          const IntVect ivFine = iv * m_refRat + IntVect(D_DECL(i, j, k));

          coarDataReg(iv, a_coarVar) += 1.0 / fineDataReg(ivFine, a_fineVar);
        }
      }
#if CH_SPACEDIM == 3
    }
#endif

    coarDataReg(iv, a_coarVar) = finePerCoar / coarDataReg(iv, a_coarVar);
  };

  // Kernel for doing irregular faces.
  auto irregularKernel = [&](const FaceIndex& face) -> void {
    a_coarData(face, a_coarVar) = 0.0;

    const FaceStencil& stencil = coarseningStencils(face, 0);
    for (int i = 0; i < stencil.size(); i++) {
      a_coarData(face, a_coarVar) += stencil.weight(i) / a_fineData(stencil.face(i), a_fineVar);
    }

    a_coarData(face, a_coarVar) = 1. / a_coarData(face, a_coarVar);
  };

  // Regular cells
  CH_START(t1);
  const Box& coarBox     = m_eblgCoFi.getDBL()[a_datInd];
  const Box  coarFaceBox = surroundingNodes(coarBox, a_dir);

  BoxLoops::loop(coarFaceBox, regularKernel);
  CH_STOP(t1);

  // Irregular faces
  CH_START(t2);
  FaceIterator& faceIt = m_irregFacesCoFi[a_datInd][a_dir];

  BoxLoops::loop(faceIt, irregularKernel);
  CH_STOP(t2);
}

void
EBCoarAve::conservativeAverage(EBFaceFAB&       a_coarData,
                               const EBFaceFAB& a_fineData,
                               const DataIndex& a_datInd,
                               const int&       a_coarVar,
                               const int&       a_fineVar,
                               const int&       a_dir) const noexcept
{
  CH_TIMERS("EBCoarAve::conservativeAverage(EBFaceFAB)");
  CH_TIMER("EBCoarAve::conservativeAverage(EBFaceFAB)::regular_cells", t1);
  CH_TIMER("EBCoarAve::conservativeAverage(EBFaceFAB)::irregular_cells", t2);

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_coarVar);
  CH_assert(a_fineData.nComp() > a_fineVar);

  FArrayBox&       coarDataReg = a_coarData.getFArrayBox();
  const FArrayBox& fineDataReg = a_fineData.getFArrayBox();

  const BaseIFFAB<FaceStencil>& coarseningStencils = m_faceConservativeStencils[a_datInd].at(a_dir);

  // Kernel for doing regular faces.
  const Real invFinePerCoar = 1.0 / std::pow(m_refRat, SpaceDim - 1);
  const int  xDoLoop        = (a_dir == 0) ? 0 : 1;
  const int  yDoLoop        = (a_dir == 1) ? 0 : 1;
#if CH_SPACEDIM == 3
  const int zDoLoop = (a_dir == 2) ? 0 : 1;
#endif

  auto regularKernel = [&](const IntVect& iv) -> void {
    coarDataReg(iv, a_coarVar) = 0.0;

#if CH_SPACEDIM == 3
    for (int k = 0; k <= (m_refRat - 1) * zDoLoop; k++) {
#endif
      for (int j = 0; j <= (m_refRat - 1) * yDoLoop; j++) {
        for (int i = 0; i <= (m_refRat - 1) * xDoLoop; i++) {
          const IntVect ivFine = iv * m_refRat + IntVect(D_DECL(i, j, k));

          coarDataReg(iv, a_coarVar) += fineDataReg(ivFine, a_fineVar);
        }
      }
#if CH_SPACEDIM == 3
    }
#endif

    coarDataReg(iv, a_coarVar) *= invFinePerCoar;
  };

  // Kernel for doing irregular faces.
  auto irregularKernel = [&](const FaceIndex& face) -> void {
    a_coarData(face, a_coarVar) = 0.0;

    const FaceStencil& stencil = coarseningStencils(face, 0);
    for (int i = 0; i < stencil.size(); i++) {
      a_coarData(face, a_coarVar) += stencil.weight(i) * a_fineData(stencil.face(i), a_fineVar);
    }
  };

  // Regular faces.
  CH_START(t1);
  const Box& coarBox     = m_eblgCoFi.getDBL()[a_datInd];
  const Box  coarFaceBox = surroundingNodes(coarBox, a_dir);

  BoxLoops::loop(coarFaceBox, regularKernel);
  CH_STOP(t1);

  // Irregular faces
  CH_START(t2);
  FaceIterator& faceIt = m_irregFacesCoFi[a_datInd][a_dir];

  BoxLoops::loop(faceIt, irregularKernel);
  CH_STOP(t2);
}

void
EBCoarAve::averageData(LevelData<BaseIVFAB<Real>>&       a_coarData,
                       const LevelData<BaseIVFAB<Real>>& a_fineData,
                       const Interval&                   a_variables,
                       const Average&                    a_average) const noexcept
{
  CH_TIMERS("EBCoarAve::averageData(LD<BaseIVFAB>)");
  CH_TIMER("EBCoarAve::averageData(LD<BaseIVFAB>)::average", t1);

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_variables.end());
  CH_assert(a_fineData.nComp() > a_variables.end());

  const DisjointBoxLayout& dblCoFi = m_eblgCoFi.getDBL();
  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();

  const DataIterator& dit  = dblFine.dataIterator();
  const int           nbox = dit.size();

  LevelData<BaseIVFAB<Real>> coFiData;
  BaseIVFactory<Real>        factCoFi(m_eblgCoFi.getEBISL(), m_irregSetsCoFi);
  coFiData.define(dblCoFi, a_variables.size(), IntVect::Zero, factCoFi);

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {
    CH_START(t1);
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      BaseIVFAB<Real>&       coarData = coFiData[din];
      const BaseIVFAB<Real>& fineData = a_fineData[din];

      // Switch between averaging methods.
      switch (a_average) {
      case Average::Arithmetic: {
        this->arithmeticAverage(coarData, fineData, din, 0, ivar);

        break;
      }
      case Average::Harmonic: {
        this->harmonicAverage(coarData, fineData, din, 0, ivar);

        break;
      }
      case Average::Conservative: {
        this->conservativeAverage(coarData, fineData, din, 0, ivar);

        break;
      }
      default: {
        MayDay::Error("EBCoarAve::averageData(LD<BaseIVFAB>) - logic bust");

        break;
      }
      }
    }
    CH_STOP(t1);

    Copier& ebCopier = m_ebCopiers[a_coarData.ghostVect()];
    if (!ebCopier.isDefined()) {
      ebCopier.define(m_eblgCoFi.getDBL(), m_eblgCoar.getDBL(), a_coarData.ghostVect());
    }

    coFiData.copyTo(Interval(0, 0), a_coarData, Interval(ivar, ivar), ebCopier);
  }
}

void
EBCoarAve::arithmeticAverage(BaseIVFAB<Real>&       a_coarData,
                             const BaseIVFAB<Real>& a_fineData,
                             const DataIndex&       a_datInd,
                             const int&             a_coarVar,
                             const int&             a_fineVar) const noexcept
{
  CH_TIME("EBCoarAve::arithmeticAverage(BaseIVFAB<Real>)");

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_coarVar);
  CH_assert(a_fineData.nComp() > a_fineVar);

  const BaseIVFAB<VoFStencil>& coarseningStencils = m_ebArithmeticStencils[a_datInd];

  auto irregularKernel = [&](const VolIndex& coarVoF) -> void {
    a_coarData(coarVoF, a_coarVar) = 0.0;

    const VoFStencil& stencil = coarseningStencils(coarVoF, 0);
    for (int i = 0; i < stencil.size(); i++) {
      a_coarData(coarVoF, a_coarVar) += stencil.weight(i) * a_fineData(stencil.vof(i), a_fineVar);
    }
  };

  VoFIterator& vofit = m_irregCellsCoFi[a_datInd];

  BoxLoops::loop(vofit, irregularKernel);
}

void
EBCoarAve::harmonicAverage(BaseIVFAB<Real>&       a_coarData,
                           const BaseIVFAB<Real>& a_fineData,
                           const DataIndex&       a_datInd,
                           const int&             a_coarVar,
                           const int&             a_fineVar) const noexcept
{
  CH_TIME("EBCoarAve::harmonicAverage(BaseIVFAB<Real>)");

  const BaseIVFAB<VoFStencil>& coarseningStencils = m_ebHarmonicStencils[a_datInd];

  auto irregularKernel = [&](const VolIndex& coarVoF) -> void {
    a_coarData(coarVoF, a_coarVar) = 0.0;

    const VoFStencil& stencil = coarseningStencils(coarVoF, 0);
    for (int i = 0; i < stencil.size(); i++) {
      a_coarData(coarVoF, a_coarVar) += stencil.weight(i) / a_fineData(stencil.vof(i), a_fineVar);
    }

    a_coarData(coarVoF, a_coarVar) = 1.0 / a_coarData(coarVoF, a_coarVar);
  };

  VoFIterator& vofit = m_irregCellsCoFi[a_datInd];

  BoxLoops::loop(vofit, irregularKernel);
}

void
EBCoarAve::conservativeAverage(BaseIVFAB<Real>&       a_coarData,
                               const BaseIVFAB<Real>& a_fineData,
                               const DataIndex&       a_datInd,
                               const int&             a_coarVar,
                               const int&             a_fineVar) const noexcept
{
  CH_TIME("EBCoarAve::conservativeAverage(BaseIVFAB<Real>)");

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_coarVar);
  CH_assert(a_fineData.nComp() > a_fineVar);

  const BaseIVFAB<VoFStencil>& coarseningStencils = m_ebConservativeStencils[a_datInd];

  auto irregularKernel = [&](const VolIndex& coarVoF) -> void {
    a_coarData(coarVoF, a_coarVar) = 0.0;

    const VoFStencil& stencil = coarseningStencils(coarVoF, 0);
    for (int i = 0; i < stencil.size(); i++) {
      a_coarData(coarVoF, a_coarVar) += stencil.weight(i) * a_fineData(stencil.vof(i), a_fineVar);
    }
  };

  VoFIterator& vofit = m_irregCellsCoFi[a_datInd];

  BoxLoops::loop(vofit, irregularKernel);
}

#include <CD_NamespaceFooter.H>
