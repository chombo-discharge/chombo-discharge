/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbCoarAve.cpp
  @brief  Implementation of CD_EbCoarAve.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBCellFactory.H>
#include <EBFluxFactory.H>
#include <BaseIVFactory.H>
#include <EBAverageF_F.H>

// Our includes
#include <CD_EbCoarAve.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

#if 1 // Debug cod eto be removed
void
EbCoarAve::averageFaceData(LevelData<EBFluxFAB>&       a_coarData,
                           const LevelData<EBFluxFAB>& a_fineData,
                           const Interval&             a_variables)
{
  this->averageData(a_coarData, a_fineData, a_variables, Average::Conservative);
}
#endif

EbCoarAve::EbCoarAve() { EBCoarseAverage::setDefaultValues(); }

EbCoarAve::~EbCoarAve() {}

EbCoarAve::EbCoarAve(const DisjointBoxLayout& a_dblFine,
                     const DisjointBoxLayout& a_dblCoar,
                     const EBISLayout&        a_ebislFine,
                     const EBISLayout&        a_ebislCoar,
                     const ProblemDomain&     a_domainCoar,
                     const int&               a_refRat,
                     const int&               a_nComp,
                     const EBIndexSpace*      a_ebisPtr)
{
  CH_TIME("EbCoarAve::EbCoarAve");

  EBCoarseAverage::setDefaultValues();
  EBCoarseAverage::define(a_dblFine, a_dblCoar, a_ebislFine, a_ebislCoar, a_domainCoar, a_refRat, a_nComp, a_ebisPtr);
}

EbCoarAve::EbCoarAve(const EBLevelGrid& a_eblgFine,
                     const EBLevelGrid& a_eblgCoar,
                     const EBLevelGrid& a_eblgCoFi,
                     const int&         a_refRat,
                     const int&         a_nComp)
{
  CH_TIME("EbCoarAve::EbCoarAve");

  EBCoarseAverage::setDefaultValues();
  EBCoarseAverage::define(a_eblgFine, a_eblgCoar, a_eblgCoFi, a_refRat, a_nComp);
}

void
EbCoarAve::averageData(LevelData<EBCellFAB>&       a_coarData,
                       const LevelData<EBCellFAB>& a_fineData,
                       const Interval&             a_variables,
                       const Average&              a_average)
{
  CH_TIME("EbCoarAve::averageData(LD<EBCellFAB>)");

  CH_assert(isDefined());

  LevelData<EBCellFAB> coarFiData;
  LevelData<EBCellFAB> fineBuffer;

  // Get input data onto the buffer if we must
  const int     ncomp = a_variables.end() + 1;
  EBCellFactory factCoFi(m_eblgCoFi.getEBISL());
  coarFiData.define(m_eblgCoFi.getDBL(), ncomp, IntVect::Zero, factCoFi);
  if (m_useFineBuffer) {
    EBCellFactory factFine(m_eblgFine.getEBISL());
    fineBuffer.define(m_eblgFine.getDBL(), ncomp, IntVect::Zero, factFine);
  }

  if (m_useFineBuffer) {
    a_fineData.copyTo(a_variables, fineBuffer, a_variables);
  }

  for (DataIterator dit = m_eblgFine.getDBL().dataIterator(); dit.ok(); ++dit) {
    const EBCellFAB* fineFAB = nullptr;
    if (m_useFineBuffer) {
      fineFAB = &fineBuffer[dit()];
    }
    else {
      fineFAB = &a_fineData[dit()];
    }

    // Switch between methods.
    switch (a_average) {
    case Average::Arithmetic: {
      this->arithmeticAverage(coarFiData[dit()], *fineFAB, dit(), a_variables, a_variables);

      break;
    }
    case Average::Harmonic: {
      this->harmonicAverage(coarFiData[dit()], *fineFAB, dit(), a_variables, a_variables);

      break;
    }
    case Average::Conservative: {
      this->conservativeAverage(coarFiData[dit()], *fineFAB, dit(), a_variables, a_variables);

      break;
    }
    default: {
      MayDay::Error("EbCoarAve::averageData(LD<EBCellFAB>) - logic bust");

      break;
    }
    }
  }

  coarFiData.copyTo(a_variables, a_coarData, a_variables);
}

void
EbCoarAve::arithmeticAverage(EBCellFAB&       a_coarData,
                             const EBCellFAB& a_fineData,
                             const DataIndex& a_datInd,
                             const Interval&  a_fineInterv,
                             const Interval&  a_coarInterv)
{
  CH_TIME("EBCoarseAverage::arithmeticAverage(EBCellFAB)");

  CH_assert(isDefined());
  CH_assert(a_fineInterv.size() == a_coarInterv.size());

  const Box& coarBox = m_eblgCoFi.getDBL()[a_datInd];
  const Box  refiBox = Box(IntVect::Zero, (m_refRat - 1) * IntVect::Unit);

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

  // Irregular cells -- recall that datInd is from the fine layout.
  const EBISBox& ebisBoxCoar = m_eblgCoFi.getEBISL()[a_datInd];
  const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[a_datInd];

  const IntVectSet& coarIrregIVS = ebisBoxCoar.getIrregIVS(coarBox);

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
EbCoarAve::harmonicAverage(EBCellFAB&       a_coarData,
                           const EBCellFAB& a_fineData,
                           const DataIndex& a_datInd,
                           const Interval&  a_fineInterv,
                           const Interval&  a_coarInterv)
{
  CH_TIME("EBCoarseAverage::arithmeticAverage(EBCellFAB)");

  CH_assert(isDefined());
  CH_assert(a_fineInterv.size() == a_coarInterv.size());

  const Box& coarBox = m_eblgCoFi.getDBL()[a_datInd];
  const Box  refiBox = Box(IntVect::Zero, (m_refRat - 1) * IntVect::Unit);

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

  // Irregular cells -- recall that datInd is from the fine layout.
  const EBISBox& ebisBoxCoar = m_eblgCoFi.getEBISL()[a_datInd];
  const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[a_datInd];

  const IntVectSet& coarIrregIVS = ebisBoxCoar.getIrregIVS(coarBox);

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
EbCoarAve::conservativeAverage(EBCellFAB&       a_coarData,
                               const EBCellFAB& a_fineData,
                               const DataIndex& a_datInd,
                               const Interval&  a_fineInterv,
                               const Interval&  a_coarInterv)
{
  CH_TIME("EBCoarseAverage::conservativeAverage(EBCellFAB)");

  CH_assert(isDefined());
  CH_assert(a_fineInterv.size() == a_coarInterv.size());

  const Real dxCoar = 1.0;
  const Real dxFine = dxCoar / m_refRat;

  const Box& coarBox = m_eblgCoFi.getDBL()[a_datInd];
  const Box  refiBox = Box(IntVect::Zero, (m_refRat - 1) * IntVect::Unit);

  // Regular cells.
  for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
    BaseFab<Real>&       coarDataReg = a_coarData.getSingleValuedFAB();
    const BaseFab<Real>& fineDataReg = a_fineData.getSingleValuedFAB();

    const int fineVar = a_fineInterv.begin() + ioff;
    const int coarVar = a_coarInterv.begin() + ioff;

    auto regularKernel = [&](const IntVect& iv) -> void {
      coarDataReg(iv, coarVar) = 0.0;

      for (BoxIterator bit(refiBox); bit.ok(); ++bit) {
        coarDataReg(iv, coarVar) += fineDataReg(m_refRat * iv + bit(), fineVar);
      }

      coarDataReg(iv, coarVar) *= std::pow(dxFine / dxCoar, SpaceDim);
    };

    BoxLoops::loop(coarBox, regularKernel);
  }

  // Irregular cells -- recall that datInd is from the fine layout.
  const EBISBox& ebisBoxCoar = m_eblgCoFi.getEBISL()[a_datInd];
  const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[a_datInd];

  const IntVectSet& coarIrregIVS = ebisBoxCoar.getIrregIVS(coarBox);

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
EbCoarAve::averageData(LevelData<EBFluxFAB>&       a_coarData,
                       const LevelData<EBFluxFAB>& a_fineData,
                       const Interval&             a_variables,
                       const Average&              a_average)
{
  CH_TIME("EbCoarAve::averageData(LD<EBFluxFAB>)");

  CH_assert(isDefined());

  LevelData<EBFluxFAB> coarFiData;
  LevelData<EBFluxFAB> fineBuffer;

  EBFluxFactory factCoFi(m_eblgCoFi.getEBISL());
  EBFluxFactory factFine(m_eblgFine.getEBISL());

  // Define buffer data.
  const int numComp = a_variables.end() + 1;
  coarFiData.define(m_eblgCoFi.getDBL(), numComp, IntVect::Zero, factCoFi);
  if (m_useFineBuffer) {
    fineBuffer.define(m_eblgFine.getDBL(), numComp, IntVect::Zero, factFine);
    a_fineData.copyTo(a_variables, fineBuffer, a_variables);
  }

  for (DataIterator dit = m_eblgFine.getDBL().dataIterator(); dit.ok(); ++dit) {
    const EBFluxFAB* fineFABPtr = nullptr;

    if (m_useFineBuffer) {
      fineFABPtr = &fineBuffer[dit()];
    }
    else {
      fineFABPtr = &a_fineData[dit()];
    }

    EBFluxFAB&       cofiFAB = coarFiData[dit()];
    const EBFluxFAB& fineFAB = *fineFABPtr;

    switch (a_average) {
    case Average::Arithmetic: {
      for (int idir = 0; idir < SpaceDim; idir++) {
        this->arithmeticAverage(cofiFAB[idir], fineFAB[idir], dit(), a_variables, a_variables, idir);
      }

      break;
    }
    case Average::Harmonic: {
      for (int idir = 0; idir < SpaceDim; idir++) {
        this->harmonicAverage(cofiFAB[idir], fineFAB[idir], dit(), a_variables, a_variables, idir);
      }

      break;
    }
    case Average::Conservative: {
      for (int idir = 0; idir < SpaceDim; idir++) {
        this->conservativeAverage(cofiFAB[idir], fineFAB[idir], dit(), a_variables, a_variables, idir);
      }

      break;
    }
    default: {
      MayDay::Error("EbCoarAve::averageData(LD<EBFluxFAB>) - logic bust");

      break;
    }
    }
  }

  // Copy back to input data.
  coarFiData.copyTo(a_variables, a_coarData, a_variables);
}

void
EbCoarAve::arithmeticAverage(EBFaceFAB&       a_coarData,
                             const EBFaceFAB& a_fineData,
                             const DataIndex& a_datInd,
                             const Interval&  a_fineInterv,
                             const Interval&  a_coarInterv,
                             const int&       a_dir)
{
  CH_TIME("EBCoarseAverage::arithmeticAverage");

  CH_assert(isDefined());
  CH_assert(a_fineInterv.size() == a_coarInterv.size());

  const int finePerCoar = std::pow(m_refRat, SpaceDim - 1);

  const Box& coarBox = m_eblgCoFi.getDBL()[a_datInd];

  Box coarFaceBox = coarBox;
  coarFaceBox.surroundingNodes(a_dir);

  const int xDoLoop = (a_dir == 0) ? 0 : 1;
  const int yDoLoop = (a_dir == 1) ? 0 : 1;
  const int zDoLoop = (a_dir == 2) ? 0 : 1;

  // Regular cells.
  for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
    const int fineVar = a_fineInterv.begin() + ioff;
    const int coarVar = a_coarInterv.begin() + ioff;

    BaseFab<Real>&       coarDataReg = a_coarData.getSingleValuedFAB();
    const BaseFab<Real>& fineDataReg = a_fineData.getSingleValuedFAB();

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

  // Now do the irregular faces.
  const EBISBox&   ebisBoxCoar  = m_eblgCoFi.getEBISL()[a_datInd];
  const EBGraph&   ebgraphCoar  = ebisBoxCoar.getEBGraph();
  const IntVectSet coarIrregIVS = ebisBoxCoar.getIrregIVS(coarBox);

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
EbCoarAve::harmonicAverage(EBFaceFAB&       a_coarData,
                           const EBFaceFAB& a_fineData,
                           const DataIndex& a_datInd,
                           const Interval&  a_fineInterv,
                           const Interval&  a_coarInterv,
                           const int&       a_dir)
{
  CH_TIME("EBCoarseAverage::harmonicAverage");

  CH_assert(isDefined());
  CH_assert(a_fineInterv.size() == a_coarInterv.size());

  const Real dxCoar      = 1.0;
  const Real dxFine      = dxCoar / m_refRat;
  const int  finePerCoar = std::pow(m_refRat, SpaceDim - 1);

  const Box& coarBox = m_eblgCoFi.getDBL()[a_datInd];

  Box coarFaceBox = coarBox;
  coarFaceBox.surroundingNodes(a_dir);

  const int xDoLoop = (a_dir == 0) ? 0 : 1;
  const int yDoLoop = (a_dir == 1) ? 0 : 1;
  const int zDoLoop = (a_dir == 2) ? 0 : 1;

  // Regular cells.
  for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
    const int fineVar = a_fineInterv.begin() + ioff;
    const int coarVar = a_coarInterv.begin() + ioff;

    BaseFab<Real>&       coarDataReg = a_coarData.getSingleValuedFAB();
    const BaseFab<Real>& fineDataReg = a_fineData.getSingleValuedFAB();

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

  // Now do the irregular faces.
  const EBISBox&   ebisBoxCoar  = m_eblgCoFi.getEBISL()[a_datInd];
  const EBGraph&   ebgraphCoar  = ebisBoxCoar.getEBGraph();
  const IntVectSet coarIrregIVS = ebisBoxCoar.getIrregIVS(coarBox);

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
EbCoarAve::conservativeAverage(EBFaceFAB&       a_coarData,
                               const EBFaceFAB& a_fineData,
                               const DataIndex& a_datInd,
                               const Interval&  a_fineInterv,
                               const Interval&  a_coarInterv,
                               const int&       a_dir)
{
  CH_TIME("EBCoarseAverage::conservativeAverage");

  CH_assert(isDefined());
  CH_assert(a_fineInterv.size() == a_coarInterv.size());

  const Real dxCoar = 1.0;
  const Real dxFine = dxCoar / m_refRat;

  const int finePerCoar = std::pow(m_refRat, SpaceDim - 1);

  const Box& coarBox = m_eblgCoFi.getDBL()[a_datInd];

  Box coarFaceBox = coarBox;
  coarFaceBox.surroundingNodes(a_dir);

  const int xDoLoop = (a_dir == 0) ? 0 : 1;
  const int yDoLoop = (a_dir == 1) ? 0 : 1;
  const int zDoLoop = (a_dir == 2) ? 0 : 1;

  // Regular cells.
  for (int ioff = 0; ioff < a_fineInterv.size(); ioff++) {
    const int fineVar = a_fineInterv.begin() + ioff;
    const int coarVar = a_coarInterv.begin() + ioff;

    BaseFab<Real>&       coarDataReg = a_coarData.getSingleValuedFAB();
    const BaseFab<Real>& fineDataReg = a_fineData.getSingleValuedFAB();

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

  // Now do the irregular faces.
  const EBISBox&   ebisBoxCoar  = m_eblgCoFi.getEBISL()[a_datInd];
  const EBISBox&   ebisBoxFine  = m_eblgFine.getEBISL()[a_datInd];
  const EBGraph&   ebgraphCoar  = ebisBoxCoar.getEBGraph();
  const IntVectSet coarIrregIVS = ebisBoxCoar.getIrregIVS(coarBox);

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
EbCoarAve::averageData(LevelData<BaseIVFAB<Real>>&       a_coarData,
                       const LevelData<BaseIVFAB<Real>>& a_fineData,
                       const Interval&                   a_variables,
                       const Average&                    a_average)
{
  CH_TIME("EbCoarAve::averageData(LD<BaseIVFAB>)");

  CH_assert(isDefined());

  LevelData<BaseIVFAB<Real>> coarFiData;
  LevelData<BaseIVFAB<Real>> fineBuffer;

  BaseIVFactory<Real> factCoFi(m_eblgCoFi.getEBISL(), m_irregSetsCoFi);
  coarFiData.define(m_eblgCoFi.getDBL(), m_nComp, IntVect::Zero, factCoFi);
  if (m_useFineBuffer) {
    BaseIVFactory<Real> factFine(m_eblgFine.getEBISL(), m_irregSetsFine);
    coarFiData.define(m_eblgFine.getDBL(), m_nComp, IntVect::Zero, factFine);
  }

  if (m_useFineBuffer) {
    a_fineData.copyTo(a_variables, fineBuffer, a_variables);
  }

  for (DataIterator dit = m_eblgFine.getDBL().dataIterator(); dit.ok(); ++dit) {
    const BaseIVFAB<Real>* fineFABPtr = nullptr;
    if (m_useFineBuffer) {
      fineFABPtr = &fineBuffer[dit()];
    }
    else {
      fineFABPtr = &a_fineData[dit()];
    }
    BaseIVFAB<Real>&       cofiFAB = coarFiData[dit()];
    const BaseIVFAB<Real>& fineFAB = *fineFABPtr;

    // Switch between averaging methods.
    switch (a_average) {
    case Average::Arithmetic: {
      this->arithmeticAverage(a_coarData[dit()], a_fineData[dit()], dit(), a_variables, a_variables);

      break;
    }
    case Average::Harmonic: {
      this->harmonicAverage(a_coarData[dit()], a_fineData[dit()], dit(), a_variables, a_variables);

      break;
    }
    case Average::Conservative: {
      this->conservativeAverage(a_coarData[dit()], a_fineData[dit()], dit(), a_variables, a_variables);

      break;
    }
    default: {
      MayDay::Error("EbCoarAve::averageData(LD<BaseIVFAB>) - logic bust");

      break;
    }
    }
  }

  coarFiData.copyTo(a_variables, a_coarData, a_variables);
}

void
EbCoarAve::arithmeticAverage(BaseIVFAB<Real>&       a_coarData,
                             const BaseIVFAB<Real>& a_fineData,
                             const DataIndex&       a_datInd,
                             const Interval&        a_coarInterv,
                             const Interval&        a_fineInterv) const
{
  CH_assert(isDefined());

  const EBISBox& ebisBoxCoar = m_eblgCoFi.getEBISL()[a_datInd];
  const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[a_datInd];

  const IntVectSet& coarIrregIVS = a_coarData.getIVS();
  const IntVectSet& fineIrregIVS = a_fineData.getIVS();

  for (VoFIterator vofitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph()); vofitCoar.ok(); ++vofitCoar) {
    const VolIndex&        coarVoF  = vofitCoar();
    const Vector<VolIndex> fineVoFs = m_eblgCoFi.getEBISL().refine(coarVoF, m_refRat, a_datInd);

    const int numFineVoFs = fineVoFs.size();

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
EbCoarAve::harmonicAverage(BaseIVFAB<Real>&       a_coarData,
                           const BaseIVFAB<Real>& a_fineData,
                           const DataIndex&       a_datInd,
                           const Interval&        a_coarInterv,
                           const Interval&        a_fineInterv) const
{
  CH_assert(isDefined());

  const EBISBox& ebisBoxCoar = m_eblgCoFi.getEBISL()[a_datInd];
  const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[a_datInd];

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
EbCoarAve::conservativeAverage(BaseIVFAB<Real>&       a_coarData,
                               const BaseIVFAB<Real>& a_fineData,
                               const DataIndex&       a_datInd,
                               const Interval&        a_coarInterv,
                               const Interval&        a_fineInterv) const
{
  CH_assert(isDefined());

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

          a_coarData(coarVoF, coarVar) += ebisBoxFine.bndryArea(fineVoF) * a_fineData(fineVoF, fineVar);
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
