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

EBCoarAve::EBCoarAve() { m_isDefined = false; }

EBCoarAve::~EBCoarAve() {}

EBCoarAve::EBCoarAve(const DisjointBoxLayout& a_dblFine,
                     const DisjointBoxLayout& a_dblCoar,
                     const EBISLayout&        a_ebislFine,
                     const EBISLayout&        a_ebislCoar,
                     const ProblemDomain&     a_domainCoar,
                     const int&               a_refRat,
                     const EBIndexSpace*      a_ebisPtr)
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
                     const int&         a_refRat)
{
  CH_TIME("EBCoarAve::EBCoarAve(EBLevelGrid version)");

  m_isDefined = false;

  this->define(a_eblgFine, a_eblgCoar, a_eblgCoFi, a_refRat);
}

void
EBCoarAve::define(const EBLevelGrid& a_eblgFine,
                  const EBLevelGrid& a_eblgCoar,
                  const EBLevelGrid& a_eblgCoFi,
                  const int&         a_refRat)
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
                       const Average&              a_average)
{
  CH_TIME("EBCoarAve::averageData(LD<EBCellFAB>)");

  CH_assert(m_isDefined);

  const Interval buffInterv = Interval(0, a_variables.size() - 1);
  const Interval fineInterv = a_variables;

  LevelData<EBCellFAB> coarFiData;
  EBCellFactory        factCoFi(m_eblgCoFi.getEBISL());
  coarFiData.define(m_eblgCoFi.getDBL(), a_variables.size(), IntVect::Zero, factCoFi);

  const DisjointBoxLayout& dbl = m_eblgFine.getDBL();

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    EBCellFAB&       coarData = coarFiData[dit()];
    const EBCellFAB& fineData = a_fineData[dit()];

    // Switch between methods.
    switch (a_average) {
    case Average::Arithmetic: {
      this->arithmeticAverage(coarData, fineData, dit(), buffInterv, fineInterv);

      break;
    }
    case Average::Harmonic: {
      this->harmonicAverage(coarData, fineData, dit(), buffInterv, fineInterv);

      break;
    }
    case Average::Conservative: {
      this->conservativeAverage(coarData, fineData, dit(), buffInterv, fineInterv);

      break;
    }
    default: {
      MayDay::Error("EBCoarAve::averageData(LD<EBCellFAB>) - logic bust");

      break;
    }
    }
  }

  coarFiData.copyTo(buffInterv, a_coarData, a_variables);
}

void
EBCoarAve::arithmeticAverage(EBCellFAB&       a_coarData,
                             const EBCellFAB& a_fineData,
                             const DataIndex& a_datInd,
                             const Interval&  a_coarInterv,
                             const Interval&  a_fineInterv)
{
  CH_TIME("EBCoarAve::arithmeticAverage(EBCellFAB)");

  CH_assert(m_isDefined);
  CH_assert(a_fineInterv.size() == a_coarInterv.size());

  const Box&        coarBox      = m_eblgCoFi.getDBL()[a_datInd];
  const Box         refiBox      = Box(IntVect::Zero, (m_refRat - 1) * IntVect::Unit);
  const EBISBox&    ebisBoxCoar  = m_eblgCoFi.getEBISL()[a_datInd];
  const IntVectSet& coarIrregIVS = ebisBoxCoar.getIrregIVS(coarBox);

  // Regular cells.
  for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
    BaseFab<Real>&       coarDataReg = a_coarData.getSingleValuedFAB();
    const BaseFab<Real>& fineDataReg = a_fineData.getSingleValuedFAB();

    const int fineVar = a_fineInterv.begin() + ioff;
    const int coarVar = a_coarInterv.begin() + ioff;

    auto regularKernel = [&](const IntVect& iv) -> void {
      Real coarVal = 0.0;

      for (BoxIterator bit(refiBox); bit.ok(); ++bit) {
        coarVal += fineDataReg(m_refRat * iv + bit(), fineVar);
      }

      coarDataReg(iv, coarVar) = coarVal / refiBox.numPts();
    };

    BoxLoops::loop(coarBox, regularKernel);
  }

  // Irregular cells.
  for (VoFIterator vofitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph()); vofitCoar.ok(); ++vofitCoar) {
    const VolIndex& coarVoF = vofitCoar();

    const Vector<VolIndex> fineVoFs = m_eblgCoFi.getEBISL().refine(coarVoF, m_refRat, a_datInd);

    const int numFineVoFs = fineVoFs.size();

    for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
      const int fineVar = a_fineInterv.begin() + ioff;
      const int coarVar = a_coarInterv.begin() + ioff;

      a_coarData(coarVoF, coarVar) = 0.0;

      // Set phic = sum(phif)/numFine
      if (numFineVoFs > 0) {
        for (int ifine = 0; ifine < numFineVoFs; ifine++) {
          const VolIndex& fineVoF = fineVoFs[ifine];

          a_coarData(coarVoF, coarVar) += a_fineData(fineVoF, fineVar);
        }

        a_coarData(coarVoF, coarVar) = a_coarData(coarVoF, coarVar) / numFineVoFs;
      }
    }
  }
}

void
EBCoarAve::harmonicAverage(EBCellFAB&       a_coarData,
                           const EBCellFAB& a_fineData,
                           const DataIndex& a_datInd,
                           const Interval&  a_coarInterv,
                           const Interval&  a_fineInterv)

{
  CH_TIME("EBCoarAve::arithmeticAverage(EBCellFAB)");

  CH_assert(m_isDefined);
  CH_assert(a_fineInterv.size() == a_coarInterv.size());

  const Box&        coarBox      = m_eblgCoFi.getDBL()[a_datInd];
  const Box         refiBox      = Box(IntVect::Zero, (m_refRat - 1) * IntVect::Unit);
  const EBISBox&    ebisBoxCoar  = m_eblgCoFi.getEBISL()[a_datInd];
  const EBISBox&    ebisBoxFine  = m_eblgFine.getEBISL()[a_datInd];
  const IntVectSet& coarIrregIVS = ebisBoxCoar.getIrregIVS(coarBox);

  // Regular cells.
  for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
    BaseFab<Real>&       coarDataReg = a_coarData.getSingleValuedFAB();
    const BaseFab<Real>& fineDataReg = a_fineData.getSingleValuedFAB();

    const int fineVar = a_fineInterv.begin() + ioff;
    const int coarVar = a_coarInterv.begin() + ioff;

    auto regularKernel = [&](const IntVect& iv) -> void {
      Real coarVal = 0.0;

      for (BoxIterator bit(refiBox); bit.ok(); ++bit) {
        coarVal += 1.0 / fineDataReg(m_refRat * iv + bit(), fineVar);
      }

      coarDataReg(iv, coarVar) *= refiBox.numPts() / coarVal;
    };

    BoxLoops::loop(coarBox, regularKernel);
  }

  // Irregular cells
  for (VoFIterator vofitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph()); vofitCoar.ok(); ++vofitCoar) {
    const VolIndex& coarVoF = vofitCoar();

    const Vector<VolIndex> fineVoFs = m_eblgCoFi.getEBISL().refine(coarVoF, m_refRat, a_datInd);

    const int numFineVoFs = fineVoFs.size();

    for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
      const int fineVar = a_fineInterv.begin() + ioff;
      const int coarVar = a_coarInterv.begin() + ioff;

      a_coarData(coarVoF, coarVar) = 0.0;

      // Set phic = numFine/sum(1/phif)
      if (numFineVoFs > 0) {
        for (int ifine = 0; ifine < numFineVoFs; ifine++) {
          const VolIndex& fineVoF = fineVoFs[ifine];

          a_coarData(coarVoF, coarVar) += 1.0 / a_fineData(fineVoF, fineVar);
        }

        a_coarData(coarVoF, coarVar) = numFineVoFs / a_coarData(coarVoF, coarVar);
      }
    }
  }
}

void
EBCoarAve::conservativeAverage(EBCellFAB&       a_coarData,
                               const EBCellFAB& a_fineData,
                               const DataIndex& a_datInd,
                               const Interval&  a_coarInterv,
                               const Interval&  a_fineInterv)
{
  CH_TIME("EBCoarAve::conservativeAverage(EBCellFAB)");

  CH_assert(m_isDefined);
  CH_assert(a_fineInterv.size() == a_coarInterv.size());

  const Real        dxCoar       = 1.0;
  const Real        dxFine       = dxCoar / m_refRat;
  const Box&        coarBox      = m_eblgCoFi.getDBL()[a_datInd];
  const Box         refiBox      = Box(IntVect::Zero, (m_refRat - 1) * IntVect::Unit);
  const EBISBox&    ebisBoxCoar  = m_eblgCoFi.getEBISL()[a_datInd];
  const EBISBox&    ebisBoxFine  = m_eblgFine.getEBISL()[a_datInd];
  const IntVectSet& coarIrregIVS = ebisBoxCoar.getIrregIVS(coarBox);

  // Regular cells.
  for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
    BaseFab<Real>&       coarDataReg = a_coarData.getSingleValuedFAB();
    const BaseFab<Real>& fineDataReg = a_fineData.getSingleValuedFAB();

    const int coarVar = a_coarInterv.begin() + ioff;
    const int fineVar = a_fineInterv.begin() + ioff;

    auto regularKernel = [&](const IntVect& iv) -> void {
      coarDataReg(iv, coarVar) = 0.0;

      for (BoxIterator bit(refiBox); bit.ok(); ++bit) {
        coarDataReg(iv, coarVar) += fineDataReg(m_refRat * iv + bit(), fineVar);
      }

      coarDataReg(iv, coarVar) *= std::pow(dxFine / dxCoar, SpaceDim);
    };

    BoxLoops::loop(coarBox, regularKernel);
  }

  // Irregular cells.

  for (VoFIterator vofitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph()); vofitCoar.ok(); ++vofitCoar) {
    const VolIndex& coarVoF = vofitCoar();
    const Real      kappaC  = ebisBoxCoar.volFrac(coarVoF);

    const Vector<VolIndex> fineVoFs = m_eblgCoFi.getEBISL().refine(coarVoF, m_refRat, a_datInd);

    const int numFineVoFs = fineVoFs.size();

    for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
      const int fineVar = a_fineInterv.begin() + ioff;
      const int coarVar = a_coarInterv.begin() + ioff;

      a_coarData(coarVoF, coarVar) = 0.0;

      // Set phic = sum(kappaf * phif * dxf^D)/(kappac * dxc^D) if we can.
      if (numFineVoFs > 0) {
        if (kappaC > 0.0) {
          for (int ifine = 0; ifine < numFineVoFs; ifine++) {
            const VolIndex& fineVoF = fineVoFs[ifine];

            a_coarData(coarVoF, coarVar) += ebisBoxFine.volFrac(fineVoF) * a_fineData(fineVoF, fineVar);
          }

          a_coarData(coarVoF, coarVar) = a_coarData(coarVoF, coarVar) * std::pow(dxFine / dxCoar, SpaceDim);
          a_coarData(coarVoF, coarVar) = a_coarData(coarVoF, coarVar) / kappaC;
        }
        else {
          // No real volume so take the arithmetic average.
          for (int ifine = 0; ifine < numFineVoFs; ifine++) {
            a_coarData(coarVoF, coarVar) += a_fineData(fineVoFs[ifine], fineVar);
          }

          a_coarData(coarVoF, coarVar) = a_coarData(coarVoF, coarVar) / numFineVoFs;
        }
      }
    }
  }
}

void
EBCoarAve::averageData(LevelData<EBFluxFAB>&       a_coarData,
                       const LevelData<EBFluxFAB>& a_fineData,
                       const Interval&             a_variables,
                       const Average&              a_average)
{
  CH_TIME("EBCoarAve::averageData(LD<EBFluxFAB>)");

  CH_assert(m_isDefined);

  const Interval buffInterv = Interval(0, a_variables.size() - 1);
  const Interval fineInterv = a_variables;

  LevelData<EBFluxFAB> coarFiData;
  EBFluxFactory        factCoFi(m_eblgCoFi.getEBISL());
  coarFiData.define(m_eblgCoFi.getDBL(), a_variables.size(), IntVect::Zero, factCoFi);

  for (DataIterator dit = m_eblgFine.getDBL().dataIterator(); dit.ok(); ++dit) {
    EBFluxFAB&       coarData = coarFiData[dit()];
    const EBFluxFAB& fineData = a_fineData[dit()];

    switch (a_average) {
    case Average::Arithmetic: {
      for (int dir = 0; dir < SpaceDim; dir++) {
        this->arithmeticAverage(coarData[dir], fineData[dir], dit(), buffInterv, fineInterv, dir);
      }

      break;
    }
    case Average::Harmonic: {
      for (int dir = 0; dir < SpaceDim; dir++) {
        this->harmonicAverage(coarData[dir], fineData[dir], dit(), buffInterv, fineInterv, dir);
      }

      break;
    }
    case Average::Conservative: {
      for (int dir = 0; dir < SpaceDim; dir++) {
        this->conservativeAverage(coarData[dir], fineData[dir], dit(), buffInterv, fineInterv, dir);
      }

      break;
    }
    default: {
      MayDay::Error("EBCoarAve::averageData(LD<EBFluxFAB>) - logic bust");

      break;
    }
    }
  }

  // Copy back to input data.
  coarFiData.copyTo(buffInterv, a_coarData, fineInterv);
}

void
EBCoarAve::arithmeticAverage(EBFaceFAB&       a_coarData,
                             const EBFaceFAB& a_fineData,
                             const DataIndex& a_datInd,
                             const Interval&  a_fineInterv,
                             const Interval&  a_coarInterv,
                             const int&       a_dir)
{
  CH_TIME("EBCoarAve::arithmeticAverage");

  CH_assert(m_isDefined);
  CH_assert(a_fineInterv.size() == a_coarInterv.size());

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
#ifdef CH_SPACEDIM == 3
  const int zDoLoop     = (a_dir == 2) ? 0 : 1;
#endif

  // Regular faces
  for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
    const int fineVar = a_fineInterv.begin() + ioff;
    const int coarVar = a_coarInterv.begin() + ioff;

    auto regularKernel = [&](const IntVect& iv) -> void {
      coarDataReg(iv, coarVar) = 0.0;

#if CH_SPACEDIM == 3
      for (int k = 0; k <= (m_refRat - 1) * zDoLoop; k++) {
#endif
        for (int j = 0; j <= (m_refRat - 1) * yDoLoop; j++) {
          for (int i = 0; i <= (m_refRat - 1) * xDoLoop; i++) {
            const IntVect ivFine = iv * m_refRat + IntVect(D_DECL(i, j, k));

            coarDataReg(iv, coarVar) += fineDataReg(ivFine, fineVar);
          }
        }
#if CH_SPACEDIM == 3
      }
#endif

      coarDataReg(iv, coarVar) *= 1.0 / finePerCoar;
    };

    BoxLoops::loop(coarFaceBox, regularKernel);
  }

  // Irregular faces
  FaceIterator faceIt(coarIrregIVS, ebgraphCoar, a_dir, FaceStop::SurroundingWithBoundary);

  for (faceIt.reset(); faceIt.ok(); ++faceIt) {
    const FaceIndex& coarFace = faceIt();

    const std::vector<FaceIndex> fineFaces = m_eblgCoFi.getEBISL().refine(coarFace, m_refRat, a_datInd).stdVector();

    const int numFineFaces = fineFaces.size();

    for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
      const int fineVar = a_fineInterv.begin() + ioff;
      const int coarVar = a_coarInterv.begin() + ioff;

      a_coarData(coarFace, coarVar) = 0.0;

      if (numFineFaces > 0) {
        for (const auto& fineFace : fineFaces) {
          a_coarData(coarFace, coarVar) += a_fineData(fineFace, fineVar);
        }

        a_coarData(coarFace, coarVar) *= 1.0 / numFineFaces;
      }
    }
  }
}

void
EBCoarAve::harmonicAverage(EBFaceFAB&       a_coarData,
                           const EBFaceFAB& a_fineData,
                           const DataIndex& a_datInd,
                           const Interval&  a_fineInterv,
                           const Interval&  a_coarInterv,
                           const int&       a_dir)
{
  CH_TIME("EBCoarAve::harmonicAverage");

  CH_assert(m_isDefined);
  CH_assert(a_fineInterv.size() == a_coarInterv.size());

  const Real dxCoar = 1.0;

  const Box&       coarBox      = m_eblgCoFi.getDBL()[a_datInd];
  const Box        coarFaceBox  = surroundingNodes(coarBox, a_dir);
  const EBISBox&   ebisBoxCoar  = m_eblgCoFi.getEBISL()[a_datInd];
  const EBGraph&   ebgraphCoar  = ebisBoxCoar.getEBGraph();
  const IntVectSet coarIrregIVS = ebisBoxCoar.getIrregIVS(coarBox);

  const int finePerCoar = std::pow(m_refRat, SpaceDim - 1);
  const int xDoLoop     = (a_dir == 0) ? 0 : 1;
  const int yDoLoop     = (a_dir == 1) ? 0 : 1;
  const int zDoLoop     = (a_dir == 2) ? 0 : 1;

  BaseFab<Real>&       coarDataReg = a_coarData.getSingleValuedFAB();
  const BaseFab<Real>& fineDataReg = a_fineData.getSingleValuedFAB();

  // Regular cells.
  for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
    const int fineVar = a_fineInterv.begin() + ioff;
    const int coarVar = a_coarInterv.begin() + ioff;

    auto regularKernel = [&](const IntVect& iv) -> void {
      coarDataReg(iv, coarVar) = 0.0;

#if CH_SPACEDIM == 3
      for (int k = 0; k <= (m_refRat - 1) * zDoLoop; k++) {
#endif
        for (int j = 0; j <= (m_refRat - 1) * yDoLoop; j++) {
          for (int i = 0; i <= (m_refRat - 1) * xDoLoop; i++) {
            const IntVect ivFine = iv * m_refRat + IntVect(D_DECL(i, j, k));

            coarDataReg(iv, coarVar) += 1.0 / fineDataReg(ivFine, fineVar);
          }
        }
#if CH_SPACEDIM == 3
      }
#endif

      coarDataReg(iv, coarVar) = finePerCoar / coarDataReg(iv, coarVar);
      ;
    };

    BoxLoops::loop(coarFaceBox, regularKernel);
  }

  // Irregular faces.
  FaceIterator faceIt(coarIrregIVS, ebgraphCoar, a_dir, FaceStop::SurroundingWithBoundary);

  for (faceIt.reset(); faceIt.ok(); ++faceIt) {
    const FaceIndex& coarFace = faceIt();

    const std::vector<FaceIndex> fineFaces = m_eblgCoFi.getEBISL().refine(coarFace, m_refRat, a_datInd).stdVector();

    const int numFineFaces = fineFaces.size();

    for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
      const int fineVar = a_fineInterv.begin() + ioff;
      const int coarVar = a_coarInterv.begin() + ioff;

      a_coarData(coarFace, coarVar) = 0.0;

      if (numFineFaces > 0) {
        for (const auto& fineFace : fineFaces) {
          a_coarData(coarFace, coarVar) += 1.0 / a_fineData(fineFace, fineVar);
        }

        a_coarData(coarFace, coarVar) = numFineFaces / a_coarData(coarFace, coarVar);
      }
    }
  }
}

void
EBCoarAve::conservativeAverage(EBFaceFAB&       a_coarData,
                               const EBFaceFAB& a_fineData,
                               const DataIndex& a_datInd,
                               const Interval&  a_fineInterv,
                               const Interval&  a_coarInterv,
                               const int&       a_dir)
{
  CH_TIME("EBCoarAve::conservativeAverage");

  CH_assert(m_isDefined);
  CH_assert(a_fineInterv.size() == a_coarInterv.size());

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
  const int zDoLoop     = (a_dir == 2) ? 0 : 1;

  BaseFab<Real>&       coarDataReg = a_coarData.getSingleValuedFAB();
  const BaseFab<Real>& fineDataReg = a_fineData.getSingleValuedFAB();

  // Regular cells.
  for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
    const int fineVar = a_fineInterv.begin() + ioff;
    const int coarVar = a_coarInterv.begin() + ioff;

    auto regularKernel = [&](const IntVect& iv) -> void {
      coarDataReg(iv, coarVar) = 0.0;

#if CH_SPACEDIM == 3
      for (int k = 0; k <= (m_refRat - 1) * zDoLoop; k++) {
#endif
        for (int j = 0; j <= (m_refRat - 1) * yDoLoop; j++) {
          for (int i = 0; i <= (m_refRat - 1) * xDoLoop; i++) {
            const IntVect ivFine = iv * m_refRat + IntVect(D_DECL(i, j, k));

            coarDataReg(iv, coarVar) += fineDataReg(ivFine, fineVar);
          }
        }
#if CH_SPACEDIM == 3
      }
#endif

      coarDataReg(iv, coarVar) = coarDataReg(iv, coarVar) / finePerCoar;
    };

    BoxLoops::loop(coarFaceBox, regularKernel);
  }

  // Irregular faces
  FaceIterator faceIt(coarIrregIVS, ebgraphCoar, a_dir, FaceStop::SurroundingWithBoundary);

  for (faceIt.reset(); faceIt.ok(); ++faceIt) {
    const FaceIndex& coarFace = faceIt();

    const std::vector<FaceIndex> fineFaces = m_eblgCoFi.getEBISL().refine(coarFace, m_refRat, a_datInd).stdVector();

    const Real areaCoar = ebisBoxCoar.areaFrac(coarFace);

    for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
      const int fineVar = a_fineInterv.begin() + ioff;
      const int coarVar = a_coarInterv.begin() + ioff;

      a_coarData(coarFace, coarVar) = 0.0;

      // Set phiCoar = sum(areaFine * phiFine) / areaCoar * (dxFine/dxCoar)^D-1
      if (areaCoar > 0.0) {
        for (const auto& fineFace : fineFaces) {
          a_coarData(coarFace, coarVar) += ebisBoxFine.areaFrac(fineFace) * a_fineData(fineFace, fineVar);
        }

        a_coarData(coarFace, coarVar) = a_coarData(coarFace, coarVar) * std::pow(dxFine / dxCoar, SpaceDim - 1);
        a_coarData(coarFace, coarVar) = a_coarData(coarFace, coarVar) / areaCoar;
      }
      else {
        const int numFineFaces = fineFaces.size();

        if (numFineFaces > 0) {
          for (const auto& fineFace : fineFaces) {
            a_coarData(coarFace, coarVar) += a_fineData(fineFace, fineVar);
          }

          a_coarData(coarFace, coarVar) = a_coarData(coarFace, coarVar) / numFineFaces;
        }
      }
    }
  }
}

void
EBCoarAve::averageData(LevelData<BaseIVFAB<Real>>&       a_coarData,
                       const LevelData<BaseIVFAB<Real>>& a_fineData,
                       const Interval&                   a_variables,
                       const Average&                    a_average)
{
  CH_TIME("EBCoarAve::averageData(LD<BaseIVFAB>)");

  CH_assert(m_isDefined);

  LevelData<BaseIVFAB<Real>> coarFiData;

  const Interval buffInterv = Interval(0, a_variables.size() - 1);
  const Interval fineInterv = a_variables;

  BaseIVFactory<Real> factCoFi(m_eblgCoFi.getEBISL(), m_irregSetsCoFi);
  coarFiData.define(m_eblgCoFi.getDBL(), a_variables.size(), IntVect::Zero, factCoFi);

  for (DataIterator dit = m_eblgFine.getDBL().dataIterator(); dit.ok(); ++dit) {
    BaseIVFAB<Real>&       coarData = coarFiData[dit()];
    const BaseIVFAB<Real>& fineData = a_fineData[dit()];

    // Switch between averaging methods.
    switch (a_average) {
    case Average::Arithmetic: {
      this->arithmeticAverage(coarData, fineData, dit(), buffInterv, fineInterv);

      break;
    }
    case Average::Harmonic: {
      this->harmonicAverage(coarData, fineData, dit(), buffInterv, fineInterv);

      break;
    }
    case Average::Conservative: {
      this->conservativeAverage(coarData, fineData, dit(), buffInterv, fineInterv);

      break;
    }
    default: {
      MayDay::Error("EBCoarAve::averageData(LD<BaseIVFAB>) - logic bust");

      break;
    }
    }
  }

  coarFiData.copyTo(buffInterv, a_coarData, a_variables);
}

void
EBCoarAve::arithmeticAverage(BaseIVFAB<Real>&       a_coarData,
                             const BaseIVFAB<Real>& a_fineData,
                             const DataIndex&       a_datInd,
                             const Interval&        a_coarInterv,
                             const Interval&        a_fineInterv) const
{
  CH_TIME("EBCoarAve::arithmeticAverage(BaseIVFAB<Real>)");

  CH_assert(m_isDefined);
  CH_assert(a_coarInterv.size() == a_fineInterv.size());

  const EBISBox& ebisBoxCoar = m_eblgCoFi.getEBISL()[a_datInd];

  const IntVectSet& coarIrregIVS = a_coarData.getIVS();
  const IntVectSet& fineIrregIVS = a_fineData.getIVS();

  for (VoFIterator vofitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph()); vofitCoar.ok(); ++vofitCoar) {
    const VolIndex&        coarVoF  = vofitCoar();
    const Vector<VolIndex> fineVoFs = m_eblgCoFi.getEBISL().refine(coarVoF, m_refRat, a_datInd);

    for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
      const int fineVar = a_fineInterv.begin() + ioff;
      const int coarVar = a_coarInterv.begin() + ioff;

      a_coarData(coarVoF, coarVar) = 0.0;

      int numVoFs = 0;

      for (int ifine = 0; ifine < fineVoFs.size(); ifine++) {
        const VolIndex& fineVoF = fineVoFs[ifine];

        if (fineIrregIVS.contains(fineVoF.gridIndex())) {
          numVoFs += 1;
          a_coarData(coarVoF, coarVar) += a_fineData(fineVoF, fineVar);
        }
      }

      if (numVoFs > 0) {
        a_coarData(coarVoF, coarVar) = a_coarData(coarVoF, coarVar) / numVoFs;
      }
    }
  }
}

void
EBCoarAve::harmonicAverage(BaseIVFAB<Real>&       a_coarData,
                           const BaseIVFAB<Real>& a_fineData,
                           const DataIndex&       a_datInd,
                           const Interval&        a_coarInterv,
                           const Interval&        a_fineInterv) const
{
  CH_TIME("EBCoarAve::harmonicAverage(BaseIVFAB<Real>)");

  CH_assert(m_isDefined);
  CH_assert(a_coarInterv.size() == a_fineInterv.size());

  const EBISBox& ebisBoxCoar = m_eblgCoFi.getEBISL()[a_datInd];

  const IntVectSet& coarIrregIVS = a_coarData.getIVS();
  const IntVectSet& fineIrregIVS = a_fineData.getIVS();

  for (VoFIterator vofitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph()); vofitCoar.ok(); ++vofitCoar) {
    const VolIndex&        coarVoF  = vofitCoar();
    const Vector<VolIndex> fineVoFs = m_eblgCoFi.getEBISL().refine(coarVoF, m_refRat, a_datInd);

    for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
      const int fineVar = a_fineInterv.begin() + ioff;
      const int coarVar = a_coarInterv.begin() + ioff;

      a_coarData(coarVoF, coarVar) = 0.0;

      int numVoFs = 0;

      for (int ifine = 0; ifine < fineVoFs.size(); ifine++) {
        const VolIndex& fineVoF = fineVoFs[ifine];

        if (fineIrregIVS.contains(fineVoF.gridIndex())) {
          numVoFs += 1;

          a_coarData(coarVoF, coarVar) += 1.0 / a_fineData(fineVoF, fineVar);
        }
      }

      if (numVoFs > 0) {
        a_coarData(coarVoF, coarVar) = numVoFs / a_coarData(coarVoF, coarVar);
      }
    }
  }
}

void
EBCoarAve::conservativeAverage(BaseIVFAB<Real>&       a_coarData,
                               const BaseIVFAB<Real>& a_fineData,
                               const DataIndex&       a_datInd,
                               const Interval&        a_coarInterv,
                               const Interval&        a_fineInterv) const
{
  CH_TIME("EBCoarAve::conservativeAverage(BaseIVFAB<Real>)");

  CH_assert(m_isDefined);
  CH_assert(a_coarInterv.size() == a_fineInterv.size());

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

    for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
      const int fineVar = a_fineInterv.begin() + ioff;
      const int coarVar = a_coarInterv.begin() + ioff;

      a_coarData(coarVoF, coarVar) = 0.0;

      if (areaCoar > 0.0) {
        for (int ifine = 0; ifine < fineVoFs.size(); ifine++) {
          const VolIndex& fineVoF = fineVoFs[ifine];

          if (fineIrregIVS.contains(fineVoF.gridIndex())) {
            a_coarData(coarVoF, coarVar) += ebisBoxFine.bndryArea(fineVoF) * a_fineData(fineVoF, fineVar);
          }
        }

        a_coarData(coarVoF, coarVar) = a_coarData(coarVoF, coarVar) * std::pow(dxFine / dxCoar, SpaceDim - 1);
        a_coarData(coarVoF, coarVar) = a_coarData(coarVoF, coarVar) / areaCoar;
      }
      else {
        // Do arithmetic average if there is no coarse boundary area.
        int numVoFs = 0;

        for (int ifine = 0; ifine < fineVoFs.size(); ifine++) {
          const VolIndex& fineVoF = fineVoFs[ifine];

          if (fineIrregIVS.contains(fineVoF.gridIndex())) {
            numVoFs += 1;

            a_coarData(coarVoF, coarVar) += a_fineData(fineVoF, fineVar);
          }
        }
        if (numVoFs > 0) {
          a_coarData(coarVoF, coarVar) = a_coarData(coarVoF, coarVar) / numVoFs;
        }
      }
    }
  }
}

#include <CD_NamespaceFooter.H>
