/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
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

EBCoarAve::~EBCoarAve() noexcept { CH_TIME("EBCoarAve::~EBCoarAve"); }

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

  m_eblgFine = a_eblgFine;
  m_eblgCoar = a_eblgCoar;
  m_eblgCoFi = a_eblgCoFi;

  m_refRat = a_refRat;

  m_irregSetsCoFi.define(m_eblgCoFi.getDBL());
  for (DataIterator dit = m_eblgCoFi.getDBL().dataIterator(); dit.ok(); ++dit) {
    m_irregSetsCoFi[dit()] = m_eblgCoFi.getEBISL()[dit()].getIrregIVS(m_eblgCoFi.getDBL().get(dit()));
  }

  m_isDefined = true;
}

void
EBCoarAve::averageData(LevelData<EBCellFAB>&       a_coarData,
                       const LevelData<EBCellFAB>& a_fineData,
                       const Interval&             a_variables,
                       const Average&              a_average) const noexcept
{
  CH_TIMERS("EBCoarAve::averageData(ebcellfab_no_buffer)");
  CH_TIMER("EBCoarAve::define_buffer", t1);
  CH_TIMER("EBCoarAve::averageData", t2);

  CH_assert(m_isDefined);
  CH_assert(a_fineData.nComp() > a_variables.end());
  CH_assert(a_coarData.nComp() > a_variables.end());

  CH_START(t1);
  LevelData<EBCellFAB> coarFiData;
  EBCellFactory        factCoFi(m_eblgCoFi.getEBISL());
  coarFiData.define(m_eblgCoFi.getDBL(), 1, IntVect::Zero, factCoFi);
  CH_STOP(t1);

  CH_START(t2);
  this->averageData(a_coarData, coarFiData, a_fineData, a_variables, a_average);
  CH_STOP(t2);
}

void
EBCoarAve::averageData(LevelData<EBCellFAB>&       a_coarData,
                       LevelData<EBCellFAB>&       a_coFiData,
                       const LevelData<EBCellFAB>& a_fineData,
                       const Interval&             a_variables,
                       const Average&              a_average) const noexcept
{
  CH_TIMERS("EBCoarAve::averageData(LD<EBCellFAB>)");
  CH_TIMER("EBCoarAve::average", t1);
  CH_TIMER("EBCoarAve::copyTo", t2);

  CH_assert(m_isDefined);
  CH_assert(a_fineData.nComp() > a_variables.end());
  CH_assert(a_coarData.nComp() > a_variables.end());

  const Interval buffInterv = Interval(0, a_variables.size() - 1);
  const Interval fineInterv = a_variables;

  const DisjointBoxLayout& dbl = m_eblgFine.getDBL();

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {
    CH_START(t1);
    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      EBCellFAB&       coarData = a_coFiData[dit()];
      const EBCellFAB& fineData = a_fineData[dit()];

      // Switch between methods.
      switch (a_average) {
      case Average::Arithmetic: {
        this->arithmeticAverage(coarData, fineData, dit(), 0, ivar);

        break;
      }
      case Average::Harmonic: {
        this->harmonicAverage(coarData, fineData, dit(), 0, ivar);

        break;
      }
      case Average::Conservative: {
        this->conservativeAverage(coarData, fineData, dit(), 0, ivar);

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
    a_coFiData.copyTo(Interval(0, 0), a_coarData, Interval(ivar, ivar));
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
  CH_TIMER("EBCoarAve::regular_cells", t1);
  CH_TIMER("EBCoarAve::irregular_cells", t2);

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_coarVar);
  CH_assert(a_fineData.nComp() > a_fineVar);

  // Regular cells.
  CH_START(t1);
  const Box&        coarBox      = m_eblgCoFi.getDBL()[a_datInd];
  const Box         refiBox      = Box(IntVect::Zero, (m_refRat - 1) * IntVect::Unit);
  const EBISBox&    ebisBoxCoar  = m_eblgCoFi.getEBISL()[a_datInd];
  const IntVectSet& coarIrregIVS = ebisBoxCoar.getIrregIVS(coarBox);

  BaseFab<Real>&       coarDataReg = a_coarData.getSingleValuedFAB();
  const BaseFab<Real>& fineDataReg = a_fineData.getSingleValuedFAB();

  auto regularKernel = [&](const IntVect& iv) -> void {
    Real coarVal = 0.0;

    for (BoxIterator bit(refiBox); bit.ok(); ++bit) {
      coarVal += fineDataReg(m_refRat * iv + bit(), a_fineVar);
    }

    coarDataReg(iv, a_coarVar) = coarVal / refiBox.numPts();
  };

  BoxLoops::loop(coarBox, regularKernel);
  CH_STOP(t1);

  // Irregular cells.
  CH_START(t2);
  for (VoFIterator vofitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph()); vofitCoar.ok(); ++vofitCoar) {
    const VolIndex& coarVoF = vofitCoar();

    const Vector<VolIndex> fineVoFs = m_eblgCoFi.getEBISL().refine(coarVoF, m_refRat, a_datInd);

    const int numFineVoFs = fineVoFs.size();

    a_coarData(coarVoF, a_coarVar) = 0.0;

    // Set phic = sum(phif)/numFine
    if (numFineVoFs > 0) {
      for (int ifine = 0; ifine < numFineVoFs; ifine++) {
        const VolIndex& fineVoF = fineVoFs[ifine];

        a_coarData(coarVoF, a_coarVar) += a_fineData(fineVoF, a_fineVar);
      }

      a_coarData(coarVoF, a_coarVar) = a_coarData(coarVoF, a_coarVar) / numFineVoFs;
    }
  }
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
  CH_TIMER("EBCoarAve::regular_cells", t1);
  CH_TIMER("EBCoarAve::irregular_cells", t2);

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_coarVar);
  CH_assert(a_fineData.nComp() > a_fineVar);

  // Regular cells
  CH_START(t1);
  const Box&        coarBox      = m_eblgCoFi.getDBL()[a_datInd];
  const Box         refiBox      = Box(IntVect::Zero, (m_refRat - 1) * IntVect::Unit);
  const EBISBox&    ebisBoxCoar  = m_eblgCoFi.getEBISL()[a_datInd];
  const IntVectSet& coarIrregIVS = ebisBoxCoar.getIrregIVS(coarBox);

  BaseFab<Real>&       coarDataReg = a_coarData.getSingleValuedFAB();
  const BaseFab<Real>& fineDataReg = a_fineData.getSingleValuedFAB();

  auto regularKernel = [&](const IntVect& iv) -> void {
    Real coarVal = 0.0;

    for (BoxIterator bit(refiBox); bit.ok(); ++bit) {
      coarVal += 1.0 / fineDataReg(m_refRat * iv + bit(), a_fineVar);
    }

    coarDataReg(iv, a_coarVar) *= refiBox.numPts() / coarVal;
  };

  BoxLoops::loop(coarBox, regularKernel);
  CH_STOP(t1);

  // Irregular cells
  CH_START(t2);
  for (VoFIterator vofitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph()); vofitCoar.ok(); ++vofitCoar) {
    const VolIndex&        coarVoF     = vofitCoar();
    const Vector<VolIndex> fineVoFs    = m_eblgCoFi.getEBISL().refine(coarVoF, m_refRat, a_datInd);
    const int              numFineVoFs = fineVoFs.size();

    a_coarData(coarVoF, a_coarVar) = 0.0;

    // Set phic = numFine/sum(1/phif)
    if (numFineVoFs > 0) {
      for (int ifine = 0; ifine < numFineVoFs; ifine++) {
        const VolIndex& fineVoF = fineVoFs[ifine];

        a_coarData(coarVoF, a_coarVar) += 1.0 / a_fineData(fineVoF, a_fineVar);
      }

      a_coarData(coarVoF, a_coarVar) = numFineVoFs / a_coarData(coarVoF, a_coarVar);
    }
  }
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
  CH_TIMER("EBCoarAve::regular_cells", t1);
  CH_TIMER("EBCoarAve::irregular_cells", t2);

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_coarVar);
  CH_assert(a_fineData.nComp() > a_fineVar);

  // Regular cells
  CH_START(t1);
  const Real        dxCoar       = 1.0;
  const Real        dxFine       = dxCoar / m_refRat;
  const Box&        coarBox      = m_eblgCoFi.getDBL()[a_datInd];
  const Box         refiBox      = Box(IntVect::Zero, (m_refRat - 1) * IntVect::Unit);
  const EBISBox&    ebisBoxCoar  = m_eblgCoFi.getEBISL()[a_datInd];
  const EBISBox&    ebisBoxFine  = m_eblgFine.getEBISL()[a_datInd];
  const IntVectSet& coarIrregIVS = ebisBoxCoar.getIrregIVS(coarBox);

  BaseFab<Real>&       coarDataReg = a_coarData.getSingleValuedFAB();
  const BaseFab<Real>& fineDataReg = a_fineData.getSingleValuedFAB();

  auto regularKernel = [&](const IntVect& iv) -> void {
    coarDataReg(iv, a_coarVar) = 0.0;

    for (BoxIterator bit(refiBox); bit.ok(); ++bit) {
      coarDataReg(iv, a_coarVar) += fineDataReg(m_refRat * iv + bit(), a_fineVar);
    }

    coarDataReg(iv, a_coarVar) *= std::pow(dxFine / dxCoar, SpaceDim);
  };

  BoxLoops::loop(coarBox, regularKernel);
  CH_STOP(t1);

  // Irregular cells.
  CH_START(t2);
  for (VoFIterator vofitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph()); vofitCoar.ok(); ++vofitCoar) {
    const VolIndex&        coarVoF     = vofitCoar();
    const Real             kappaC      = ebisBoxCoar.volFrac(coarVoF);
    const Vector<VolIndex> fineVoFs    = m_eblgCoFi.getEBISL().refine(coarVoF, m_refRat, a_datInd);
    const int              numFineVoFs = fineVoFs.size();

    a_coarData(coarVoF, a_coarVar) = 0.0;

    // Set phic = sum(kappaf * phif * dxf^D)/(kappac * dxc^D) if we can.
    if (numFineVoFs > 0) {
      if (kappaC > 0.0) {
        for (int ifine = 0; ifine < numFineVoFs; ifine++) {
          const VolIndex& fineVoF = fineVoFs[ifine];

          a_coarData(coarVoF, a_coarVar) += ebisBoxFine.volFrac(fineVoF) * a_fineData(fineVoF, a_fineVar);
        }

        a_coarData(coarVoF, a_coarVar) = a_coarData(coarVoF, a_coarVar) * std::pow(dxFine / dxCoar, SpaceDim);
        a_coarData(coarVoF, a_coarVar) = a_coarData(coarVoF, a_coarVar) / kappaC;
      }
      else {
        // No real volume so take the arithmetic average.
        for (int ifine = 0; ifine < numFineVoFs; ifine++) {
          a_coarData(coarVoF, a_coarVar) += a_fineData(fineVoFs[ifine], a_fineVar);
        }

        a_coarData(coarVoF, a_coarVar) = a_coarData(coarVoF, a_coarVar) / numFineVoFs;
      }
    }
  }
  CH_STOP(t2);
}

void
EBCoarAve::averageData(LevelData<EBFluxFAB>&       a_coarData,
                       const LevelData<EBFluxFAB>& a_fineData,
                       const Interval&             a_variables,
                       const Average&              a_average) const noexcept
{
  CH_TIMERS("EBCoarAve::averageData(ebfluxfab_no_buffer)");
  CH_TIMER("EBCoarAve::define_buffer", t1);
  CH_TIMER("EBCoarAve::averageData", t2);

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_variables.end());
  CH_assert(a_fineData.nComp() > a_variables.end());

  CH_START(t1);
  LevelData<EBFluxFAB> coarFiData;
  EBFluxFactory        factCoFi(m_eblgCoFi.getEBISL());
  coarFiData.define(m_eblgCoFi.getDBL(), 1, IntVect::Zero, factCoFi);
  CH_STOP(t1);

  CH_START(t2);
  this->averageData(a_coarData, coarFiData, a_fineData, a_variables, a_average);
  CH_STOP(t2);
}

void
EBCoarAve::averageData(LevelData<EBFluxFAB>&       a_coarData,
                       LevelData<EBFluxFAB>&       a_coFiData,
                       const LevelData<EBFluxFAB>& a_fineData,
                       const Interval&             a_variables,
                       const Average&              a_average) const noexcept
{
  CH_TIMERS("EBCoarAve::averageData(LD<EBFluxFAB>)");
  CH_TIMER("EBCoarAve::average", t1);
  CH_TIMER("EBCoarAve::copyTo", t2);

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_variables.end());
  CH_assert(a_fineData.nComp() > a_variables.end());

  CH_START(t1);
  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {
    for (DataIterator dit = m_eblgFine.getDBL().dataIterator(); dit.ok(); ++dit) {
      EBFluxFAB&       coarData = a_coFiData[dit()];
      const EBFluxFAB& fineData = a_fineData[dit()];

      switch (a_average) {
      case Average::Arithmetic: {
        for (int dir = 0; dir < SpaceDim; dir++) {
          this->arithmeticAverage(coarData[dir], fineData[dir], dit(), 0, ivar, dir);
        }

        break;
      }
      case Average::Harmonic: {
        for (int dir = 0; dir < SpaceDim; dir++) {
          this->harmonicAverage(coarData[dir], fineData[dir], dit(), 0, ivar, dir);
        }

        break;
      }
      case Average::Conservative: {
        for (int dir = 0; dir < SpaceDim; dir++) {
          this->conservativeAverage(coarData[dir], fineData[dir], dit(), 0, ivar, dir);
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

    // Copy back to input data.
    CH_START(t2);
    a_coFiData.copyTo(Interval(0, 0), a_coarData, Interval(ivar, ivar));
    CH_STOP(t2);
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
  CH_TIMER("EBCoarAve::regular_cells", t1);
  CH_TIMER("EBCoarAve::irregular_cells", t2);

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_coarVar);
  CH_assert(a_fineData.nComp() > a_fineVar);

  // Regular faces.
  CH_START(t1);
  const Box&       coarBox      = m_eblgCoFi.getDBL()[a_datInd];
  const Box        coarFaceBox  = surroundingNodes(coarBox, a_dir);
  const EBISBox&   ebisBoxCoar  = m_eblgCoFi.getEBISL()[a_datInd];
  const EBGraph&   ebgraphCoar  = ebisBoxCoar.getEBGraph();
  const IntVectSet coarIrregIVS = ebisBoxCoar.getIrregIVS(coarBox);

  BaseFab<Real>&       coarDataReg = a_coarData.getSingleValuedFAB();
  const BaseFab<Real>& fineDataReg = a_fineData.getSingleValuedFAB();

  const int finePerCoar = std::pow(m_refRat, SpaceDim - 1);
  const int xDoLoop     = (a_dir == 0) ? 0 : 1;
  const int yDoLoop     = (a_dir == 1) ? 0 : 1;
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

    coarDataReg(iv, a_coarVar) *= 1.0 / finePerCoar;
  };

  BoxLoops::loop(coarFaceBox, regularKernel);
  CH_STOP(t1);

  // Irregular faces
  CH_START(t2);
  FaceIterator faceIt(coarIrregIVS, ebgraphCoar, a_dir, FaceStop::SurroundingWithBoundary);

  for (faceIt.reset(); faceIt.ok(); ++faceIt) {
    const FaceIndex&             coarFace     = faceIt();
    const std::vector<FaceIndex> fineFaces    = m_eblgCoFi.getEBISL().refine(coarFace, m_refRat, a_datInd).stdVector();
    const int                    numFineFaces = fineFaces.size();

    a_coarData(coarFace, a_coarVar) = 0.0;

    if (numFineFaces > 0) {
      for (const auto& fineFace : fineFaces) {
        a_coarData(coarFace, a_coarVar) += a_fineData(fineFace, a_fineVar);
      }

      a_coarData(coarFace, a_coarVar) *= 1.0 / numFineFaces;
    }
  }
  CH_START(t2);
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
  CH_TIMER("EBCoarAve::regular_cells", t1);
  CH_TIMER("EBCoarAve::irregular_cells", t2);

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_coarVar);
  CH_assert(a_fineData.nComp() > a_fineVar);

  // Regular faces.
  CH_START(t1);
  const Box&       coarBox      = m_eblgCoFi.getDBL()[a_datInd];
  const Box        coarFaceBox  = surroundingNodes(coarBox, a_dir);
  const EBISBox&   ebisBoxCoar  = m_eblgCoFi.getEBISL()[a_datInd];
  const EBGraph&   ebgraphCoar  = ebisBoxCoar.getEBGraph();
  const IntVectSet coarIrregIVS = ebisBoxCoar.getIrregIVS(coarBox);

  const int finePerCoar = std::pow(m_refRat, SpaceDim - 1);
  const int xDoLoop     = (a_dir == 0) ? 0 : 1;
  const int yDoLoop     = (a_dir == 1) ? 0 : 1;
#if CH_SPACEDIM == 3
  const int zDoLoop = (a_dir == 2) ? 0 : 1;
#endif

  BaseFab<Real>&       coarDataReg = a_coarData.getSingleValuedFAB();
  const BaseFab<Real>& fineDataReg = a_fineData.getSingleValuedFAB();

  // Regular cells.
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

  BoxLoops::loop(coarFaceBox, regularKernel);
  CH_STOP(t1);

  // Irregular faces.
  CH_START(t2);
  FaceIterator faceIt(coarIrregIVS, ebgraphCoar, a_dir, FaceStop::SurroundingWithBoundary);

  for (faceIt.reset(); faceIt.ok(); ++faceIt) {
    const FaceIndex&             coarFace     = faceIt();
    const std::vector<FaceIndex> fineFaces    = m_eblgCoFi.getEBISL().refine(coarFace, m_refRat, a_datInd).stdVector();
    const int                    numFineFaces = fineFaces.size();

    a_coarData(coarFace, a_coarVar) = 0.0;

    if (numFineFaces > 0) {
      for (const auto& fineFace : fineFaces) {
        a_coarData(coarFace, a_coarVar) += 1.0 / a_fineData(fineFace, a_fineVar);
      }

      a_coarData(coarFace, a_coarVar) = numFineFaces / a_coarData(coarFace, a_coarVar);
    }
  }
  CH_START(t2);
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
  CH_TIMER("EBCoarAve::regular_cells", t1);
  CH_TIMER("EBCoarAve::irregular_cells", t2);

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_coarVar);
  CH_assert(a_fineData.nComp() > a_fineVar);

  CH_START(t1);
  const Real dxCoar = 1.0;
  const Real dxFine = dxCoar / m_refRat;

  const Box&       coarBox      = m_eblgCoFi.getDBL()[a_datInd];
  const Box        coarFaceBox  = surroundingNodes(coarBox, a_dir);
  const EBISBox&   ebisBoxCoar  = m_eblgCoFi.getEBISL()[a_datInd];
  const EBISBox&   ebisBoxFine  = m_eblgFine.getEBISL()[a_datInd];
  const EBGraph&   ebgraphCoar  = ebisBoxCoar.getEBGraph();
  const IntVectSet coarIrregIVS = ebisBoxCoar.getIrregIVS(coarBox);

  const int finePerCoar = std::pow(m_refRat, SpaceDim - 1);
  const int xDoLoop     = (a_dir == 0) ? 0 : 1;
  const int yDoLoop     = (a_dir == 1) ? 0 : 1;
#if CH_SPACEDIM == 3
  const int zDoLoop = (a_dir == 2) ? 0 : 1;
#endif

  BaseFab<Real>&       coarDataReg = a_coarData.getSingleValuedFAB();
  const BaseFab<Real>& fineDataReg = a_fineData.getSingleValuedFAB();

  // Regular cells.
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

    coarDataReg(iv, a_coarVar) = coarDataReg(iv, a_coarVar) / finePerCoar;
  };

  BoxLoops::loop(coarFaceBox, regularKernel);
  CH_STOP(t1);

  // Irregular faces
  CH_START(t2);
  FaceIterator faceIt(coarIrregIVS, ebgraphCoar, a_dir, FaceStop::SurroundingWithBoundary);

  for (faceIt.reset(); faceIt.ok(); ++faceIt) {
    const FaceIndex&             coarFace  = faceIt();
    const std::vector<FaceIndex> fineFaces = m_eblgCoFi.getEBISL().refine(coarFace, m_refRat, a_datInd).stdVector();
    const Real                   areaCoar  = ebisBoxCoar.areaFrac(coarFace);

    a_coarData(coarFace, a_coarVar) = 0.0;

    // Set phiCoar = sum(areaFine * phiFine) / areaCoar * (dxFine/dxCoar)^D-1
    if (areaCoar > 0.0) {
      for (const auto& fineFace : fineFaces) {
        a_coarData(coarFace, a_coarVar) += ebisBoxFine.areaFrac(fineFace) * a_fineData(fineFace, a_fineVar);
      }

      a_coarData(coarFace, a_coarVar) = a_coarData(coarFace, a_coarVar) * std::pow(dxFine / dxCoar, SpaceDim - 1);
      a_coarData(coarFace, a_coarVar) = a_coarData(coarFace, a_coarVar) / areaCoar;
    }
    else {
      const int numFineFaces = fineFaces.size();

      if (numFineFaces > 0) {
        for (const auto& fineFace : fineFaces) {
          a_coarData(coarFace, a_coarVar) += a_fineData(fineFace, a_fineVar);
        }

        a_coarData(coarFace, a_coarVar) = a_coarData(coarFace, a_coarVar) / numFineFaces;
      }
    }
  }
  CH_STOP(t2);
}

void
EBCoarAve::averageData(LevelData<BaseIVFAB<Real>>&       a_coarData,
                       const LevelData<BaseIVFAB<Real>>& a_fineData,
                       const Interval&                   a_variables,
                       const Average&                    a_average) const noexcept
{
  CH_TIMERS("EBCoarAve::averageData(baseivfab_no_buffer)");
  CH_TIMER("EBCoarAve::define_buffer", t1);
  CH_TIMER("EBCoarAve::averageData", t2);

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_variables.end());
  CH_assert(a_fineData.nComp() > a_variables.end());

  CH_START(t1);
  LevelData<BaseIVFAB<Real>> coarFiData;
  BaseIVFactory<Real>        factCoFi(m_eblgCoFi.getEBISL(), m_irregSetsCoFi);
  coarFiData.define(m_eblgCoFi.getDBL(), a_variables.size(), IntVect::Zero, factCoFi);
  CH_STOP(t1);

  CH_START(t2);
  this->averageData(a_coarData, coarFiData, a_fineData, a_variables, a_average);
  CH_STOP(t2);
}

void
EBCoarAve::averageData(LevelData<BaseIVFAB<Real>>&       a_coarData,
                       LevelData<BaseIVFAB<Real>>&       a_coFiData,
                       const LevelData<BaseIVFAB<Real>>& a_fineData,
                       const Interval&                   a_variables,
                       const Average&                    a_average) const noexcept
{
  CH_TIMERS("EBCoarAve::averageData(LD<BaseIVFAB>)");
  CH_TIMER("EBCoarAve::average", t1);
  CH_TIMER("EBCoarAve::copyTo", t2);

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_variables.end());
  CH_assert(a_fineData.nComp() > a_variables.end());

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {
    CH_START(t1);
    for (DataIterator dit = m_eblgFine.getDBL().dataIterator(); dit.ok(); ++dit) {
      BaseIVFAB<Real>&       coarData = a_coFiData[dit()];
      const BaseIVFAB<Real>& fineData = a_fineData[dit()];

      // Switch between averaging methods.
      switch (a_average) {
      case Average::Arithmetic: {
        this->arithmeticAverage(coarData, fineData, dit(), 0, ivar);

        break;
      }
      case Average::Harmonic: {
        this->harmonicAverage(coarData, fineData, dit(), 0, ivar);

        break;
      }
      case Average::Conservative: {
        this->conservativeAverage(coarData, fineData, dit(), 0, ivar);

        break;
      }
      default: {
        MayDay::Error("EBCoarAve::averageData(LD<BaseIVFAB>) - logic bust");

        break;
      }
      }
    }
    CH_STOP(t1);

    CH_START(t2);
    a_coFiData.copyTo(Interval(0, 0), a_coarData, Interval(ivar, ivar));
    CH_STOP(t2);
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

  const EBISBox& ebisBoxCoar = m_eblgCoFi.getEBISL()[a_datInd];

  const IntVectSet& coarIrregIVS = a_coarData.getIVS();
  const IntVectSet& fineIrregIVS = a_fineData.getIVS();

  for (VoFIterator vofitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph()); vofitCoar.ok(); ++vofitCoar) {
    const VolIndex&        coarVoF  = vofitCoar();
    const Vector<VolIndex> fineVoFs = m_eblgCoFi.getEBISL().refine(coarVoF, m_refRat, a_datInd);

    a_coarData(coarVoF, a_coarVar) = 0.0;

    int numVoFs = 0;

    for (int ifine = 0; ifine < fineVoFs.size(); ifine++) {
      const VolIndex& fineVoF = fineVoFs[ifine];

      if (fineIrregIVS.contains(fineVoF.gridIndex())) {
        numVoFs += 1;
        a_coarData(coarVoF, a_coarVar) += a_fineData(fineVoF, a_fineVar);
      }
    }

    if (numVoFs > 0) {
      a_coarData(coarVoF, a_coarVar) = a_coarData(coarVoF, a_coarVar) / numVoFs;
    }
  }
}

void
EBCoarAve::harmonicAverage(BaseIVFAB<Real>&       a_coarData,
                           const BaseIVFAB<Real>& a_fineData,
                           const DataIndex&       a_datInd,
                           const int&             a_coarVar,
                           const int&             a_fineVar) const noexcept
{
  CH_TIME("EBCoarAve::harmonicAverage(BaseIVFAB<Real>)");

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() > a_coarVar);
  CH_assert(a_fineData.nComp() > a_fineVar);

  const EBISBox& ebisBoxCoar = m_eblgCoFi.getEBISL()[a_datInd];

  const IntVectSet& coarIrregIVS = a_coarData.getIVS();
  const IntVectSet& fineIrregIVS = a_fineData.getIVS();

  for (VoFIterator vofitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph()); vofitCoar.ok(); ++vofitCoar) {
    const VolIndex&        coarVoF  = vofitCoar();
    const Vector<VolIndex> fineVoFs = m_eblgCoFi.getEBISL().refine(coarVoF, m_refRat, a_datInd);

    a_coarData(coarVoF, a_coarVar) = 0.0;

    int numVoFs = 0;

    for (int ifine = 0; ifine < fineVoFs.size(); ifine++) {
      const VolIndex& fineVoF = fineVoFs[ifine];

      if (fineIrregIVS.contains(fineVoF.gridIndex())) {
        numVoFs += 1;

        a_coarData(coarVoF, a_coarVar) += 1.0 / a_fineData(fineVoF, a_fineVar);
      }
    }

    if (numVoFs > 0) {
      a_coarData(coarVoF, a_coarVar) = numVoFs / a_coarData(coarVoF, a_coarVar);
    }
  }
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

  const EBISBox& ebisBoxCoar = m_eblgCoFi.getEBISL()[a_datInd];
  const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[a_datInd];

  const IntVectSet& coarIrregIVS = a_coarData.getIVS();
  const IntVectSet& fineIrregIVS = a_fineData.getIVS();

  const Real dxCoar = 1.0;
  const Real dxFine = dxCoar / m_refRat;

  // This loop computes phiCoar = sum(phiFine*areaFine)/areaCoar. If areaCoar == 0
  // then we default to an arithmetic average.
  for (VoFIterator vofitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph()); vofitCoar.ok(); ++vofitCoar) {
    const VolIndex&        coarVoF  = vofitCoar();
    const Vector<VolIndex> fineVoFs = m_eblgCoFi.getEBISL().refine(coarVoF, m_refRat, a_datInd);

    const Real areaCoar = ebisBoxCoar.bndryArea(coarVoF);

    a_coarData(coarVoF, a_coarVar) = 0.0;

    if (areaCoar > 0.0) {
      for (int ifine = 0; ifine < fineVoFs.size(); ifine++) {
        const VolIndex& fineVoF = fineVoFs[ifine];

        if (fineIrregIVS.contains(fineVoF.gridIndex())) {
          a_coarData(coarVoF, a_coarVar) += ebisBoxFine.bndryArea(fineVoF) * a_fineData(fineVoF, a_fineVar);
        }
      }

      a_coarData(coarVoF, a_coarVar) = a_coarData(coarVoF, a_coarVar) * std::pow(dxFine / dxCoar, SpaceDim - 1);
      a_coarData(coarVoF, a_coarVar) = a_coarData(coarVoF, a_coarVar) / areaCoar;
    }
    else {
      // Do arithmetic average if there is no coarse boundary area.
      int numVoFs = 0;

      for (int ifine = 0; ifine < fineVoFs.size(); ifine++) {
        const VolIndex& fineVoF = fineVoFs[ifine];

        if (fineIrregIVS.contains(fineVoF.gridIndex())) {
          numVoFs += 1;

          a_coarData(coarVoF, a_coarVar) += a_fineData(fineVoF, a_fineVar);
        }
      }
      if (numVoFs > 0) {
        a_coarData(coarVoF, a_coarVar) = a_coarData(coarVoF, a_coarVar) / numVoFs;
      }
    }
  }
}

#include <CD_NamespaceFooter.H>
