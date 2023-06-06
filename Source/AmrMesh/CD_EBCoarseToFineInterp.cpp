/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBCoarseToFineInterp.cpp
  @brief  Implementation of EBCoarseToFineInterp.H
  @author Robert Marskar
*/

// Std includes
#include <limits>

// Chombo includes
#include <CH_Timer.H>
#include <EBCellFactory.H>
#include <BaseIVFactory.H>

// Our includes
#include <CD_EBCoarseToFineInterp.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

EBCoarseToFineInterp::EBCoarseToFineInterp() noexcept
{
  CH_TIME("EBCoarseToFineInterp::EBCoarseToFineInterp(weak)");

  m_isDefined = false;
}

EBCoarseToFineInterp::EBCoarseToFineInterp(const EBLevelGrid& a_eblgFine,
                                           const EBLevelGrid& a_eblgCoFi,
                                           const EBLevelGrid& a_eblgCoar,
                                           const int&         a_refRat) noexcept
{
  CH_TIME("EBCoarseToFineInterp::EBCoarseToFineInterp(full)");

  this->define(a_eblgFine, a_eblgCoFi, a_eblgCoar, a_refRat);
}

EBCoarseToFineInterp::~EBCoarseToFineInterp() noexcept {}

void
EBCoarseToFineInterp::define(const EBLevelGrid& a_eblgFine,
                             const EBLevelGrid& a_eblgCoFi,
                             const EBLevelGrid& a_eblgCoar,
                             const int&         a_refRat) noexcept
{
  CH_TIMERS("EBCoarseToFineInterp::define");
  CH_TIMER("EBCoarseToFineInterp::define_ebpwlfineinterp", t0);
  CH_TIMER("EBCoarseToFineInterp::define_buffers", t1);

  CH_assert(a_refRat > 1);

  m_refRat   = a_refRat;
  m_eblgFine = a_eblgFine;
  m_eblgCoar = a_eblgCoar;
  m_eblgCoFi = a_eblgCoFi;

  CH_START(t1);
  // Define the irreg data holder.
  LayoutData<IntVectSet> irregFine(m_eblgFine.getDBL());
  LayoutData<IntVectSet> irregCoar(m_eblgCoFi.getDBL());

  m_fineVoFs.define(m_eblgFine.getDBL());
  m_coarVoFs.define(m_eblgFine.getDBL());

  for (DataIterator dit(m_eblgFine.getDBL()); dit.ok(); ++dit) {
    const Box& fineBox = m_eblgFine.getDBL()[dit()];
    const Box& coarBox = m_eblgCoFi.getDBL()[dit()];

    const EBISBox& fineEBISBox = m_eblgFine.getEBISL()[dit()];
    const EBISBox& coarEBISBox = m_eblgCoFi.getEBISL()[dit()];

    irregFine[dit()] = fineEBISBox.getIrregIVS(fineBox);
    irregCoar[dit()] = coarEBISBox.getIrregIVS(coarBox);

    m_fineVoFs[dit()].define(irregFine[dit()], fineEBISBox.getEBGraph());
    m_coarVoFs[dit()].define(coarEBISBox.getIrregIVS(coarBox), coarEBISBox.getEBGraph());
  }

  m_irregCoFi.define(m_eblgCoFi.getDBL(), 1, IntVect::Zero, BaseIVFactory<Real>(m_eblgCoFi.getEBISL(), irregCoar));
  m_areaWeights.define(m_eblgFine.getDBL(), 1, IntVect::Zero, BaseIVFactory<Real>(m_eblgFine.getEBISL(), irregFine));
  m_volumeWeights.define(m_eblgFine.getDBL(), 1, IntVect::Zero, BaseIVFactory<Real>(m_eblgFine.getEBISL(), irregFine));
  CH_STOP(t1);

  this->defineWeights();
  this->defineBuffers();

  m_isDefined = true;
}

void
EBCoarseToFineInterp::defineWeights() noexcept
{
  CH_TIME("EBCoarseToFineInterp::defineWeights");

  const EBISLayout& ebislFine = m_eblgFine.getEBISL();
  const EBISLayout& ebislCoar = m_eblgCoFi.getEBISL();

  const Real areaFactor   = std::pow(m_refRat, SpaceDim - 1);
  const Real volumeFactor = std::pow(m_refRat, SpaceDim);

  for (DataIterator dit(m_eblgFine.getDBL()); dit.ok(); ++dit) {
    const EBISBox& ebisBoxFine = ebislFine[dit()];
    const EBISBox& ebisBoxCoar = ebislCoar[dit()];

    BaseIVFAB<Real>& volumeWeights = m_volumeWeights[dit()];
    BaseIVFAB<Real>& areaWeights   = m_areaWeights[dit()];

    volumeWeights.setVal(0.0);
    areaWeights.setVal(0.0);

    auto kernel = [&](const VolIndex& fineVoF) -> void {
      // TLDR: We want refRat^(D-1) * coarArea/fineArea. We refine the coarse vof (which also gives regular cells)
      // but that won't matter because they don't have a boundary area. We do it similarly for the volume fractions.

      const VolIndex         coarVoF  = ebislFine.coarsen(fineVoF, m_refRat, dit());
      const Vector<VolIndex> fineVoFs = ebislCoar.refine(coarVoF, m_refRat, dit());

      Real fineVolume = 0.0;
      Real fineArea   = 0.0;

      const Real coarVolume = ebisBoxCoar.volFrac(coarVoF);
      const Real coarArea   = ebisBoxCoar.bndryArea(coarVoF);

      for (int i = 0; i < fineVoFs.size(); i++) {
        fineArea += ebisBoxFine.bndryArea(fineVoFs[i]);
        fineVolume += ebisBoxFine.volFrac(fineVoFs[i]);
      }

      if (fineArea > 0.0) {
        areaWeights(fineVoF, 0) = areaFactor * coarArea / fineArea;
      }

      if (fineVolume > 0.0) {
        volumeWeights(fineVoF, 0) = volumeFactor * coarVolume / fineVolume;
      }
    };

    BoxLoops::loop(m_fineVoFs[dit()], kernel);
  }
}

void
EBCoarseToFineInterp::defineBuffers() noexcept
{
  CH_TIME("EBCoarseToFineInterp::defineBuffers");

  // Note --same number of ghost cells as in EBCoarseToFineInterp. Must be one because of the slopes!
  m_copier.define(m_eblgCoar.getDBL(), m_eblgCoFi.getDBL(), IntVect::Unit);
}

void
EBCoarseToFineInterp::interpolate(LevelData<EBCellFAB>&             a_fineData,
                                  const LevelData<EBCellFAB>&       a_coarData,
                                  const Interval&                   a_variables,
                                  const EBCoarseToFineInterp::Type& a_interpType) const noexcept
{
  CH_TIMERS("EBCoarseToFineInterp::interpolate(LD<EBCellFAB>)");
  CH_TIMER("EBCoarseToFineInterp::interpolate(LD<EBCellFAB>)::interpolate", t1);

  CH_assert(m_isDefined);
  CH_assert(a_fineData.nComp() > a_variables.end());
  CH_assert(a_coarData.nComp() > a_variables.end());

  LevelData<EBCellFAB> m_coFiData(m_eblgCoFi.getDBL(), 1, IntVect::Unit, EBCellFactory(m_eblgCoFi.getEBISL()));

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {
    const Interval srcInterv = Interval(ivar, ivar);
    const Interval dstInterv = Interval(0, 0);

    a_coarData.copyTo(srcInterv, m_coFiData, dstInterv, m_copier);

    CH_START(t1);
    for (DataIterator dit(m_eblgCoFi.getDBL()); dit.ok(); ++dit) {
      switch (a_interpType) {
      case EBCoarseToFineInterp::Type::PWC: {
        this->interpolatePWC(a_fineData[dit()], m_coFiData[dit()], dit(), ivar, 0);

        break;
      }
      case EBCoarseToFineInterp::Type::ConservativePWC: {
        this->interpolateConservativePWC(a_fineData[dit()], m_coFiData[dit()], dit(), ivar, 0);

        break;
      }
      case EBCoarseToFineInterp::Type::ConservativeMinMod: {
        this->interpolateConservativeSlope(a_fineData[dit()], m_coFiData[dit()], dit(), ivar, 0, SlopeLimiter::MinMod);

        break;
      }
      case EBCoarseToFineInterp::Type::ConservativeMonotonizedCentral: {
        this->interpolateConservativeSlope(a_fineData[dit()],
                                           m_coFiData[dit()],
                                           dit(),
                                           ivar,
                                           0,
                                           SlopeLimiter::MonotonizedCentral);

        break;
      }
      default: {
        MayDay::Error("EBCoarseToFineInterp::interpolate - logic bust. Interpolation type not supported");

        break;
      }
      }
    }
    CH_STOP(t1);
  }
}

void
EBCoarseToFineInterp::interpolate(LevelData<BaseIVFAB<Real>>&       a_fineData,
                                  const LevelData<BaseIVFAB<Real>>& a_coarData,
                                  const Interval&                   a_variables,
                                  const EBCoarseToFineInterp::Type& a_interpType) const noexcept
{
  CH_TIMERS("EBCoarseToFineInterp::interpolate(LD<BaseIVFAB<Real>>)");
  CH_TIMER("EBCoarseToFineInterp::define_buffer", t1);
  CH_TIMER("EBCoarseToFineInterp::copyTo", t2);

  CH_assert(m_isDefined);
  CH_assert(a_fineData.nComp() > a_variables.end());
  CH_assert(a_coarData.nComp() > a_variables.end());

  const DisjointBoxLayout& dblCoFi = m_eblgCoFi.getDBL();

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {
    CH_START(t1);
    a_coarData.copyTo(Interval(ivar, ivar), m_irregCoFi, Interval(0, 0));
    CH_STOP(t1);

    for (DataIterator dit(dblCoFi); dit.ok(); ++dit) {
      switch (a_interpType) {
      case EBCoarseToFineInterp::Type::PWC: {
        this->interpolatePWC(a_fineData[dit()], m_irregCoFi[dit()], dit(), ivar, 0);

        break;
      }
      case EBCoarseToFineInterp::Type::ConservativePWC: {
        this->interpolateConservativePWC(a_fineData[dit()], m_irregCoFi[dit()], dit(), ivar, 0);

        break;
      }
      default: {
        MayDay::Error(
          "EBCoarseToFineInterp::interpolate(LD<BaseIVFAB<Real>) - logic bust. Interpolation type not supported");

        break;
      }
      }
    }
  }
}

void
EBCoarseToFineInterp::interpolatePWC(EBCellFAB&       a_fineData,
                                     const EBCellFAB& a_coarData,
                                     const DataIndex& a_dit,
                                     const int&       a_fineVar,
                                     const int&       a_coarVar) const noexcept
{
  CH_TIMERS("EBCoarseToFineInterp::interpolatePWC(EBCellFAB)");
  CH_TIMER("EBCoarseToFineInterp::regular_regrid", t1);
  CH_TIMER("EBCoarseToFineInterp::irregular_regrid", t2);

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

    for (int i = 0; i < fineVoFs.size(); i++) {
      a_fineData(fineVoFs[i], a_fineVar) = a_coarData(coarVoF, a_coarVar);
    }
  };

  CH_START(t1);
  BoxLoops::loop(coarBox, regularKernel);
  CH_STOP(t1);

  CH_START(t2);
  BoxLoops::loop(m_coarVoFs[a_dit], irregularKernel);
  CH_STOP(t2);
}

void
EBCoarseToFineInterp::interpolateConservativePWC(EBCellFAB&       a_fineData,
                                                 const EBCellFAB& a_coarData,
                                                 const DataIndex& a_dit,
                                                 const int&       a_fineVar,
                                                 const int&       a_coarVar) const noexcept
{
  CH_TIMERS("EBCoarseToFineInterp::interpolateConservativePWC(EBCellFAB)");
  CH_TIMER("EBCoarseToFineInterp::regular_regrid", t1);
  CH_TIMER("EBCoarseToFineInterp::irregular_regrid", t2);

  CH_assert(a_fineData.nComp() > a_fineVar);
  CH_assert(a_coarData.nComp() > a_coarVar);

  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();
  const DisjointBoxLayout& dblCoar = m_eblgCoFi.getDBL();

  const EBISLayout& ebislFine = m_eblgFine.getEBISL();
  const EBISLayout& ebislCoar = m_eblgCoFi.getEBISL();

  const EBISBox& ebisBoxFine = ebislFine[a_dit];
  const EBISBox& ebisBoxCoar = ebislCoar[a_dit];

  const Box fineBox = dblFine[a_dit];
  const Box coarBox = dblCoar[a_dit];

  const BaseIVFAB<Real>& volumeFactor = m_volumeWeights[a_dit];

  FArrayBox&       fineDataReg = a_fineData.getFArrayBox();
  const FArrayBox& coarDataReg = a_coarData.getFArrayBox();

  // Regular kernel. Set the fine data equal to the coarse data.
  auto regularKernel = [&](const IntVect& fineIV) -> void {
    const IntVect coarIV = coarsen(fineIV, m_refRat);

    fineDataReg(fineIV, a_fineVar) = coarDataReg(coarIV, a_coarVar);
  };

  // Cut-cell kernel.
  auto irregularKernel = [&](const VolIndex& fineVoF) -> void {
    const VolIndex& coarVoF = ebislFine.coarsen(fineVoF, m_refRat, a_dit);

    a_fineData(fineVoF, a_fineVar) = volumeFactor(fineVoF, 0) * a_coarData(coarVoF, a_coarVar);
  };

  CH_START(t1);
  BoxLoops::loop(fineBox, regularKernel);
  CH_STOP(t1);

  CH_START(t2);
  BoxLoops::loop(m_fineVoFs[a_dit], irregularKernel);
  CH_STOP(t2);

#ifndef NDEBUG
  this->checkConservation(a_fineData, a_coarData, a_dit, a_fineVar, a_coarVar);
#endif
}

void
EBCoarseToFineInterp::interpolateConservativeSlope(EBCellFAB&          a_fineData,
                                                   const EBCellFAB&    a_coarData,
                                                   const DataIndex&    a_dit,
                                                   const int&          a_fineVar,
                                                   const int&          a_coarVar,
                                                   const SlopeLimiter& a_limiter) const noexcept
{
  CH_TIMERS("EBCoarseToFineInterp::interpolateConservativeSlope(EBCellFAB)");
  CH_TIMER("EBCoarseToFineInterp::slopes_interior", t1);
  CH_TIMER("EBCoarseToFineInterp::extrap_regular", t2);
  CH_TIMER("EBCoarseToFineInterp::constant_term_reg", t3);
  CH_TIMER("EBCoarseToFineInterp::irreg_interp", t4);

  CH_assert(a_fineData.nComp() > a_fineVar);
  CH_assert(a_coarData.nComp() > a_coarVar);

  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();
  const DisjointBoxLayout& dblCoar = m_eblgCoFi.getDBL();

  const EBISLayout& ebislFine = m_eblgFine.getEBISL();
  const EBISLayout& ebislCoar = m_eblgCoFi.getEBISL();

  const ProblemDomain& domainFine = m_eblgFine.getDomain();
  const ProblemDomain& domainCoar = m_eblgCoar.getDomain();

  const EBISBox& ebisBoxFine = ebislFine[a_dit];
  const EBISBox& ebisBoxCoar = ebislCoar[a_dit];

  const Box coarBox = dblCoar[a_dit];
  const Box fineBox = dblFine[a_dit];

  const BaseIVFAB<Real>& volumeFactor = m_volumeWeights[a_dit];

  FArrayBox&       fineDataReg = a_fineData.getFArrayBox();
  const FArrayBox& coarDataReg = a_coarData.getFArrayBox();

  // Set to zero because we increment below.
  a_fineData.setVal(0.0);

  // Compute the slopes.
  EBCellFAB  slopes(ebisBoxCoar, coarBox, 1);
  FArrayBox& slopesReg = slopes.getFArrayBox();
  slopes.setVal(0.0);

  // Various slope definitions.
  auto minmod = [](const Real& dwl, const Real& dwr) -> Real {
    Real slope = 0.0;

    if (dwl * dwr > 0.0) {
      slope = std::abs(dwl) < std::abs(dwr) ? dwl : dwr;
    }
    return slope;
  };

  auto mc = [](const Real& dwl, const Real& dwr) -> Real {
    Real slope = 0.0;

    if (dwl * dwr > 0.0) {
      const Real dwc = dwl + dwr;
      const Real sgn = Real((dwc > 0.0) - (dwc < 0.0));

      slope = sgn * std::min(0.5 * std::abs(dwc), 2.0 * std::min(std::abs(dwl), std::abs(dwr)));
    }

    return slope;
  };

  auto superbee = [=](const Real& dwl, const Real& dwr) -> Real {
    Real slope = 0.0;

    if (dwl * dwr > 0.0) {
      const Real s1 = minmod(dwl, 2 * dwr);
      const Real s2 = minmod(dwr, 2 * dwl);

      if (s1 * s2 > 0.0) {
        slope = std::abs(s1) > std::abs(s2) ? s1 : s2;
      }
    }

    return slope;
  };

  // Compute slopes in every direction. Slopes in boundary cells are zero (because what would they be?). And in
  // cut-cells we set slopes to zero (because I don't know how to combine hard conservation with slope extrapolation!).
  for (int dir = 0; dir < SpaceDim; dir++) {
    const IntVect shift = BASISV(dir);

    const Box interiorCells = grow(grow(coarBox, 1) & domainCoar, -1);

    // Define the regular grid kernel.
    std::function<void(const IntVect& iv)> interiorKernel;
    switch (a_limiter) {
    case SlopeLimiter::MinMod: {
      interiorKernel = [&](const IntVect& iv) -> void {
        const Real dwl = coarDataReg(iv, a_coarVar) - coarDataReg(iv - shift, a_coarVar);
        const Real dwr = coarDataReg(iv + shift, a_coarVar) - coarDataReg(iv, a_coarVar);

        slopesReg(iv, 0) = minmod(dwl, dwr);
      };

      break;
    }
    case SlopeLimiter::MonotonizedCentral: {
      interiorKernel = [&](const IntVect& iv) -> void {
        const Real dwl = coarDataReg(iv, a_coarVar) - coarDataReg(iv - shift, a_coarVar);
        const Real dwr = coarDataReg(iv + shift, a_coarVar) - coarDataReg(iv, a_coarVar);

        slopesReg(iv, 0) = mc(dwl, dwr);
      };

      break;
    }
    case SlopeLimiter::Superbee: {
      interiorKernel = [&](const IntVect& iv) -> void {
        const Real dwl = coarDataReg(iv, a_coarVar) - coarDataReg(iv - shift, a_coarVar);
        const Real dwr = coarDataReg(iv + shift, a_coarVar) - coarDataReg(iv, a_coarVar);

        slopesReg(iv, 0) = superbee(dwl, dwr);
      };

      break;
    }
    default: {
      MayDay::Error("EBCoarseToFineInterp::interpolateConservativeSlope - logic bust 1");

      break;
    }
    }

    // Reset slopes in cut-cells because we can't extrapolate inside of them.
    auto resetSlopeIrreg = [&](const VolIndex& coarVoF) -> void {
      slopes(coarVoF, 0) = 0.0;
    };

    // Apply slopes in the fine-grid interior cells.
    auto slopeExtrapRegular = [&](const IntVect& fineIV) -> void {
      const IntVect  coarIV = coarsen(fineIV, m_refRat);
      const Real&    slope  = slopesReg(coarIV, 0);
      const RealVect delta  = (RealVect(fineIV) - m_refRat * RealVect(coarIV) + 0.5 * (1.0 - m_refRat)) / m_refRat;

      fineDataReg(fineIV, a_fineVar) += slope * delta[dir];
    };

    // Compute slopes in interior cells. Crap on boundary and cut-cells.
    CH_START(t1);
    BoxLoops::loop(interiorCells, interiorKernel);
    BoxLoops::loop(m_coarVoFs[a_dit], resetSlopeIrreg);
    CH_STOP(t1);

    CH_START(t2);
    BoxLoops::loop(fineBox, slopeExtrapRegular);
    CH_STOP(t2);
  }

  // Add in the constant term.
  auto regularConstantTerm = [&](const IntVect& fineIV) {
    const IntVect coarIV = coarsen(fineIV, m_refRat);

    fineDataReg(fineIV, a_fineVar) += coarDataReg(coarIV, a_coarVar);
  };

  // Reset irregular cells, using hard conservation.
  auto irregularInterp = [&](const VolIndex& fineVoF) {
    const VolIndex coarVoF = ebislFine.coarsen(fineVoF, m_refRat, a_dit);

    a_fineData(fineVoF, a_fineVar) = volumeFactor(fineVoF, 0) * a_coarData(coarVoF, a_coarVar);
  };

  CH_START(t3);
  BoxLoops::loop(fineBox, regularConstantTerm);
  CH_STOP(t3);

  CH_START(t4);
  BoxLoops::loop(m_fineVoFs[a_dit], irregularInterp);
  CH_STOP(t4);

#ifndef NDEBUG
  this->checkConservation(a_fineData, a_coarData, a_dit, a_fineVar, a_coarVar);
#endif
}

void
EBCoarseToFineInterp::interpolatePWC(BaseIVFAB<Real>&       a_fineData,
                                     const BaseIVFAB<Real>& a_coarData,
                                     const DataIndex&       a_dit,
                                     const int&             a_fineVar,
                                     const int&             a_coarVar) const noexcept
{
  CH_TIME("EBCoarseToFineInterp::interpolatePWC(BaseIVFAB<Real>)");

  CH_assert(a_fineData.nComp() > a_fineVar);
  CH_assert(a_coarData.nComp() > a_coarVar);

  const EBISLayout& ebislFine = m_eblgFine.getEBISL();
  const EBISLayout& ebislCoar = m_eblgCoFi.getEBISL();

  const EBISBox& fineEBISBox = ebislFine[a_dit];
  const EBISBox& coarEBISBox = ebislCoar[a_dit];

  auto kernel = [&](const VolIndex& fineVoF) -> void {
    const VolIndex coarVoF = ebislFine.coarsen(fineVoF, m_refRat, a_dit);

    CH_assert(coarEBISBox.isIrregular(coarVoF.gridIndex()));
    CH_assert(fineEBISBox.isIrregular(fineVoF.gridIndex()));

    a_fineData(fineVoF, a_fineVar) = a_coarData(coarVoF, a_coarVar);
  };

  BoxLoops::loop(m_fineVoFs[a_dit], kernel);
}

void
EBCoarseToFineInterp::interpolateConservativePWC(BaseIVFAB<Real>&       a_fineData,
                                                 const BaseIVFAB<Real>& a_coarData,
                                                 const DataIndex&       a_dit,
                                                 const int&             a_fineVar,
                                                 const int&             a_coarVar) const noexcept
{
  CH_TIME("EBCoarseToFineInterp::interpolateConservativePWC(BaseIVFAB<Real>)");

  CH_assert(a_fineData.nComp() > a_fineVar);
  CH_assert(a_coarData.nComp() > a_coarVar);

  const EBISLayout& ebislFine = m_eblgFine.getEBISL();
  const EBISLayout& ebislCoar = m_eblgCoFi.getEBISL();

  const EBISBox& fineEBISBox = ebislFine[a_dit];
  const EBISBox& coarEBISBox = ebislCoar[a_dit];

  const BaseIVFAB<Real>& weights = m_areaWeights[a_dit];

  auto kernel = [&](const VolIndex& fineVoF) -> void {
    const VolIndex coarVoF = ebislFine.coarsen(fineVoF, m_refRat, a_dit);

    CH_assert(coarEBISBox.isIrregular(coarVoF.gridIndex()));
    CH_assert(fineEBISBox.isIrregular(fineVoF.gridIndex()));

    a_fineData(fineVoF, a_fineVar) = a_coarData(coarVoF, a_coarVar) * weights(fineVoF, 0);
  };

  BoxLoops::loop(m_fineVoFs[a_dit], kernel);

#ifndef NDEBUG
  this->checkConservation(a_fineData, a_coarData, a_dit, a_fineVar, a_coarVar);
#endif
}

void
EBCoarseToFineInterp::checkConservation(const EBCellFAB& a_fineData,
                                        const EBCellFAB& a_coarData,
                                        const DataIndex& a_dit,
                                        const int        a_fineVar,
                                        const int        a_coarVar) const noexcept
{
  CH_TIME("EBCoarseToFineInterp::checkConservation(EBCellFAB)");

  CH_assert(a_fineData.nComp() > a_fineVar);
  CH_assert(a_coarData.nComp() > a_coarVar);

  const Box fineBox = m_eblgFine.getDBL()[a_dit];
  const Box coarBox = m_eblgCoFi.getDBL()[a_dit];

  const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[a_dit];
  const EBISBox& ebisBoxCoar = m_eblgCoFi.getEBISL()[a_dit];

  const FArrayBox& fineDataReg = a_fineData.getFArrayBox();
  const FArrayBox& coarDataReg = a_coarData.getFArrayBox();

  Real sumCoar = 0.0;
  Real sumFine = 0.0;

  auto regCoar = [&](const IntVect& iv) -> void {
    if (ebisBoxCoar.isRegular(iv)) {
      sumCoar += coarDataReg(iv, a_coarVar);
    }
  };
  auto regFine = [&](const IntVect& iv) -> void {
    if (ebisBoxFine.isRegular(iv)) {
      sumFine += fineDataReg(iv, a_fineVar);
    }
  };
  auto irregCoar = [&](const VolIndex& coarVoF) -> void {
    sumCoar += ebisBoxCoar.volFrac(coarVoF) * a_coarData(coarVoF, a_coarVar);
  };
  auto irregFine = [&](const VolIndex& fineVoF) -> void {
    sumFine += ebisBoxFine.volFrac(fineVoF) * a_fineData(fineVoF, a_fineVar);
  };

  BoxLoops::loop(coarBox, regCoar);
  BoxLoops::loop(fineBox, regFine);
  BoxLoops::loop(m_coarVoFs[a_dit], irregCoar);
  BoxLoops::loop(m_fineVoFs[a_dit], irregFine);

  sumCoar *= std::pow(m_refRat, SpaceDim);

  const Real delta = std::abs(sumCoar - sumFine);
  const Real eps   = std::numeric_limits<Real>::epsilon();

  // Issue an error message if we did not conserve. Try to guard against false negatives by some more
  // thorough evaluations of round-off errors.
  if (delta > eps && std::max(std::abs(sumCoar), std::abs(sumFine)) > eps &&
      delta / std::abs(std::max(sumCoar, sumFine)) > 1.E-10) {

    const std::string baseErr = "EBCoarseToFineInterp::checkConservation(EBCellFAB) - did not conserve = ";
    std::cout << baseErr << sumFine << "\t" << sumCoar << std::endl;
  }
}

void
EBCoarseToFineInterp::checkConservation(const BaseIVFAB<Real>& a_fineData,
                                        const BaseIVFAB<Real>& a_coarData,
                                        const DataIndex&       a_dit,
                                        const int              a_fineVar,
                                        const int              a_coarVar) const noexcept
{
  CH_TIME("EBCoarseToFineInterp::checkConservation(EBCellFAB)");

  CH_assert(a_fineData.nComp() > a_fineVar);
  CH_assert(a_coarData.nComp() > a_coarVar);

  const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[a_dit];
  const EBISBox& ebisBoxCoar = m_eblgCoFi.getEBISL()[a_dit];

  Real sumCoar = 0.0;
  Real sumFine = 0.0;

  auto irregCoar = [&](const VolIndex& coarVoF) -> void {
    sumCoar += ebisBoxCoar.bndryArea(coarVoF) * a_coarData(coarVoF, a_coarVar);
  };
  auto irregFine = [&](const VolIndex& fineVoF) -> void {
    sumFine += ebisBoxFine.bndryArea(fineVoF) * a_fineData(fineVoF, a_fineVar);
  };

  BoxLoops::loop(m_coarVoFs[a_dit], irregCoar);
  BoxLoops::loop(m_fineVoFs[a_dit], irregFine);

  sumCoar *= std::pow(m_refRat, SpaceDim - 1);

  const Real delta = std::abs(sumCoar - sumFine);
  const Real eps   = std::numeric_limits<Real>::epsilon();

  // Issue an error message if we did not conserve. Try to guard against false negatives by some more
  // thorough evaluations of round-off errors.
  if (delta > eps && std::max(std::abs(sumCoar), std::abs(sumFine)) > eps &&
      delta / std::abs(std::max(sumCoar, sumFine)) > 1.E-10) {

    const std::string baseErr = "EBCoarseToFineInterp::checkConservation(BaseIVFAB) - did not conserve = ";
    std::cout << baseErr << sumFine << "\t" << sumCoar << std::endl;
  }
}

#include <CD_NamespaceFooter.H>
