/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBReflux.cpp
  @brief  Implementation of CD_EBReflux.H
  @author Robert Marskar
*/

// Chombo includes
#include <NeighborIterator.H>
#include <EBCellFactory.H>
#include <EBFluxFactory.H>
#include <CH_Timer.H>

// Our includes
#include <CD_EBReflux.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

EBReflux::EBReflux() noexcept
{
  CH_TIME("EBReflux::EBReflux(weak)");

  m_isDefined = false;
}

EBReflux::EBReflux(const EBLevelGrid& a_eblg,
                   const EBLevelGrid& a_eblgFine,
                   const EBLevelGrid& a_eblgCoFi,
                   const int          a_refRat) noexcept
{
  CH_TIME("EBReflux::EBReflux(full)");

  this->define(a_eblg, a_eblgFine, a_eblgCoFi, a_refRat);
}

EBReflux::~EBReflux() noexcept { CH_TIME("EBReflux::~EBReflux"); }

void
EBReflux::define(const EBLevelGrid& a_eblg,
                 const EBLevelGrid& a_eblgFine,
                 const EBLevelGrid& a_eblgCoFi,
                 const int          a_refRat) noexcept
{
  CH_TIME("EBReflux::define");

  m_eblg     = a_eblg;
  m_eblgFine = a_eblgFine;
  m_eblgCoFi = a_eblgCoFi;
  m_refRat   = a_refRat;

  this->defineRegionsCF();
  this->defineStencils();
  this->defineBuffers();

  m_isDefined = true;
}

void
EBReflux::defineRegionsCF() noexcept
{
  CH_TIMERS("EBReflux::defineRegionsCF");
  CH_TIMER("EBReflux::defineRegionsCF::layout_define", t1);
  CH_TIMER("EBReflux::defineRegionsCF::regular_regions", t2);
  CH_TIMER("EBReflux::defineRegionsCF::irregular_regions", t3);

  const DisjointBoxLayout& dbl     = m_eblg.getDBL();
  const DisjointBoxLayout& dblCoFi = m_eblgCoFi.getDBL();

  const ProblemDomain& domain     = m_eblg.getDomain();
  const ProblemDomain& domainCoFi = m_eblgCoFi.getDomain();

  const EBISLayout& ebisl     = m_eblg.getEBISL();
  const EBISLayout& ebislCoFi = m_eblgCoFi.getEBISL();

  CH_START(t1);
  m_regularCoarseFineRegions.define(dbl);
  m_irregularCoarseFineRegions.define(dbl);

  Copier copier;
  copier.ghostDefine(dblCoFi, dbl, domain, IntVect::Unit);
  CH_STOP(t1);

  // Define the irregular coarse-fine interface. This is way more involved since we need to figure out which cells
  // on the coarse-grid side interfaces with cells on the fine-grid side.
  CH_START(t2);
  LevelData<FArrayBox> coarMask(dbl, 1, IntVect::Zero);
  LevelData<FArrayBox> coFiMask(dblCoFi, 1, IntVect::Unit);

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {

      const DataIterator& dit     = dbl.dataIterator();
      const DataIterator& ditCoFi = dblCoFi.dataIterator();

      const int nbox     = dit.size();
      const int nboxCoFi = ditCoFi.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        coarMask[din].setVal(0.0);
      }

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nboxCoFi; mybox++) {
        const DataIndex& din = ditCoFi[mybox];

        coFiMask[din].setVal(0.0);

        const Box boxCoFi     = dblCoFi[din];
        const Box sideBoxCoFi = adjCellBox(boxCoFi, dir, sit(), 1);

        if (domain.contains(sideBoxCoFi)) {
          DenseIntVectSet cfivs = DenseIntVectSet(sideBoxCoFi, true);

          NeighborIterator nit(dblCoFi);
          for (nit.begin(din); nit.ok(); ++nit) {
            cfivs -= dblCoFi[nit()];
          }

          for (DenseIntVectSetIterator ivsIt(cfivs); ivsIt.ok(); ++ivsIt) {
            coFiMask[din](ivsIt(), 0) = 1.0;
          }
        }
      }

      // Copy the data to the coarse grid.
      const Interval interv(0, 0);
      coFiMask.copyTo(interv, coarMask, interv, copier, LDaddOp<FArrayBox>());

      // Define regular cells
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        const EBISBox&   ebisBox = ebisl[din];
        const FArrayBox& mask    = coarMask[din];
        DenseIntVectSet  cfivs(dbl[din], false);

        for (BoxIterator bit(dbl[din]); bit.ok(); ++bit) {
          const IntVect iv = bit();
          if (mask(iv, 0) > 0.0 && ebisBox.isRegular(iv)) {
            cfivs |= iv;
          }
        }

        cfivs.recalcMinBox();
        m_regularCoarseFineRegions[din][std::make_pair(dir, sit())] = cfivs;
      }

      // Define irregular cells.
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        const Box        cellBox = dbl[din];
        const FArrayBox& mask    = coarMask[din];
        const EBISBox&   ebisBox = ebisl[din];
        const EBGraph&   ebGraph = ebisBox.getEBGraph();

        IntVectSet irregCells;

        auto findIrregCells = [&](const IntVect& iv) -> void {
          if (mask(iv, 0) > 0.0 && ebisBox.isIrregular(iv)) {
            irregCells |= iv;
          }
        };

        BoxLoops::loop(cellBox, findIrregCells);

        // Define appropriate iterators.
        auto& irregularCoarseFineRegions = m_irregularCoarseFineRegions[din];
        irregularCoarseFineRegions[std::make_pair(dir, sit())].define(irregCells, ebGraph);
      }
    }
  }
  CH_STOP(t2);
}

void
EBReflux::defineStencils() noexcept
{
  CH_TIME("EBReflux::defineStencils");

  const DisjointBoxLayout& dblCoar = m_eblgCoFi.getDBL();
  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();

  const DataIterator ditCoar = dblCoar.dataIterator();
  const DataIterator ditFine = dblFine.dataIterator();

  const EBISLayout& ebislCoar = m_eblgCoFi.getEBISL();
  const EBISLayout& ebislFine = m_eblgFine.getEBISL();

  const Real dxCoar   = 1.0;
  const Real dxFine   = dxCoar / m_refRat;
  const Real dxFactor = std::pow(dxFine / dxCoar, SpaceDim - 1);

  m_fluxCoarseningStencils.define(dblCoar);
  m_fluxCoarseningRegions.define(dblCoar);

  const int nbox = ditCoar.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = ditCoar[mybox];

    const Box boxCoar = dblCoar[din];
    const Box boxFine = dblFine[din];

    const EBISBox& ebisBoxCoar = ebislCoar[din];
    const EBISBox& ebisBoxFine = ebislFine[din];

    const EBGraph& graphCoar = ebisBoxCoar.getEBGraph();
    const EBGraph& graphFine = ebisBoxFine.getEBGraph();

    for (int dir = 0; dir < SpaceDim; dir++) {
      for (SideIterator sit; sit.ok(); ++sit) {
        const Box        boxSide  = adjCellBox(boxCoar, dir, sit(), -1);
        const IntVectSet irregIVS = ebisBoxCoar.getIrregIVS(boxSide);

        BaseIFFAB<FaceStencil>& faceStencils = m_fluxCoarseningStencils[din][std::make_pair(dir, sit())];
        FaceIterator&           faceIterator = m_fluxCoarseningRegions[din][std::make_pair(dir, sit())];

        faceStencils.define(irregIVS, graphCoar, dir, 1);
        faceIterator.define(irregIVS, graphCoar, dir, FaceStop::SurroundingWithBoundary);

        // Go through the coarse face and build the requires stencils for conservative flux coarsening.
        for (faceIterator.reset(); faceIterator.ok(); ++faceIterator) {
          const FaceIndex& coarFace = faceIterator();
          const Real       areaCoar = ebisBoxCoar.areaFrac(coarFace);

          FaceStencil& stencil = faceStencils(coarFace, 0);
          stencil.clear();

          if (areaCoar > 0.0) {

            const Vector<FaceIndex>& fineFaces = ebislCoar.refine(coarFace, m_refRat, din);
            for (int iface = 0; iface < fineFaces.size(); iface++) {
              const FaceIndex& fineFace   = fineFaces[iface];
              const Real       fineWeight = ebisBoxFine.areaFrac(fineFace) * dxFactor / areaCoar;

              stencil.add(fineFace, fineWeight);
            }
          }
        }
      }
    }
  }
}

void
EBReflux::defineBuffers() noexcept
{
  CH_TIME("EBReflux::defineBuffers");

  // Note: MUST have the same number of ghost cells as the buffers being defined.
  m_copier.define(m_eblgCoFi.getDBL(), m_eblg.getDBL(), IntVect::Unit);
}

void
EBReflux::reflux(LevelData<EBCellFAB>&       a_Lphi,
                 const LevelData<EBFluxFAB>& a_flux,
                 const LevelData<EBFluxFAB>& a_fineFlux,
                 const Interval              a_variables,
                 const Real                  a_scaleCoarFlux,
                 const Real                  a_scaleFineFlux) const noexcept
{
  CH_TIMERS("EBReflux::reflux");
  CH_TIMER("EBReflux::reflux::define_buffers", t1);

  CH_assert(m_isDefined);
  CH_assert(a_Lphi.nComp() > a_variables.end());
  CH_assert(a_flux.nComp() > a_variables.end());
  CH_assert(a_fineFlux.nComp() > a_variables.end());

  const DisjointBoxLayout& dbl     = m_eblg.getDBL();
  const DisjointBoxLayout& dblCoFi = m_eblgCoFi.getDBL();

  const EBISLayout& ebisl     = m_eblg.getEBISL();
  const EBISLayout& ebislCoFi = m_eblgCoFi.getEBISL();

  LevelData<EBFluxFAB> fluxCoFi(dblCoFi, 1, IntVect::Zero, EBFluxFactory(ebislCoFi));
  LevelData<EBFluxFAB> fluxCoar(dbl, 1, IntVect::Unit, EBFluxFactory(ebisl));

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {

    // Coarsen fluxes
    this->coarsenFluxesCF(fluxCoFi, a_fineFlux, 0, ivar);

    // Copy fluxes to coarse grids
    const Interval srcInterv = Interval(0, 0);
    const Interval dstInterv = Interval(0, 0);

    fluxCoFi.copyTo(srcInterv, fluxCoar, dstInterv, m_copier);

    // Reflux the coarse level.
    this->refluxIntoCoarse(a_Lphi, a_flux, fluxCoar, ivar, 0, ivar, a_scaleCoarFlux, a_scaleFineFlux);
  }
}

void
EBReflux::coarsenFluxesCF(LevelData<EBFluxFAB>&       a_coarFluxes,
                          const LevelData<EBFluxFAB>& a_fineFluxes,
                          const int                   a_coarVar,
                          const int                   a_fineVar) const noexcept
{
  CH_TIMERS("EBReflux::coarsenFluxes");
  CH_TIMER("EBReflux::coarsenFluxes::define_iterator", t1);
  CH_TIMER("EBReflux::coarsenFluxes::regular_faces", t2);
  CH_TIMER("EBReflux::coarsenFluxes::irregular_faces", t3);

  CH_assert(m_isDefined);
  CH_assert(a_coarVar >= 0);
  CH_assert(a_fineVar >= 0);
  CH_assert(a_coarFluxes.nComp() > a_coarVar);
  CH_assert(a_fineFluxes.nComp() > a_fineVar);

  const DisjointBoxLayout& dblCoar = m_eblgCoFi.getDBL();
  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();

  const DataIterator& ditCoar = dblCoar.dataIterator();

  const EBISLayout& ebislCoar = m_eblgCoFi.getEBISL();
  const EBISLayout& ebislFine = m_eblgFine.getEBISL();

  const Real dxCoar         = 1.0;
  const Real dxFine         = dxCoar / m_refRat;
  const Real invFinePerCoar = 1.0 / std::pow(m_refRat, SpaceDim - 1);

  const int nbox = ditCoar.size();
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = ditCoar[mybox];

    const Box& coarCellBox = dblCoar[din];
    const Box& fineCellBox = dblFine[din];

    const EBISBox& coarEBISBox = ebislCoar[din];
    const EBISBox& fineEBISBox = ebislFine[din];

    const EBGraph& coarGraph = coarEBISBox.getEBGraph();
    const EBGraph& fineGraph = fineEBISBox.getEBGraph();

    for (int dir = 0; dir < SpaceDim; dir++) {
      const Box coarFaceBox = surroundingNodes(coarCellBox, dir);

      EBFaceFAB&       coarFlux = a_coarFluxes[din][dir];
      const EBFaceFAB& fineFlux = a_fineFluxes[din][dir];

      FArrayBox&       coarFluxReg = coarFlux.getFArrayBox();
      const FArrayBox& fineFluxReg = fineFlux.getFArrayBox();

      // Some hooks for deciding how to do the averaging.
      const int xDoLoop = (dir == 0) ? 0 : 1;
      const int yDoLoop = (dir == 1) ? 0 : 1;
#if CH_SPACEDIM == 3
      const int zDoLoop = (dir == 2) ? 0 : 1;
#endif

      for (SideIterator sit; sit.ok(); ++sit) {

        // This is a BaseIFFAB<FaceStencil>
        const auto& coarseningStencils = m_fluxCoarseningStencils[din].at(std::make_pair(dir, sit()));

        // Kernel for regular cells.
        auto regularKernel = [&](const IntVect& iv) -> void {
          coarFluxReg(iv, a_coarVar) = 0.0;

#if CH_SPACEDIM == 3
          for (int k = 0; k <= (m_refRat - 1) * zDoLoop; k++) {
#endif
            for (int j = 0; j <= (m_refRat - 1) * yDoLoop; j++) {
              for (int i = 0; i <= (m_refRat - 1) * xDoLoop; i++) {
                const IntVect ivFine = iv * m_refRat + IntVect(D_DECL(i, j, k));

                coarFluxReg(iv, a_coarVar) += fineFluxReg(ivFine, a_fineVar);
              }
            }
#if CH_SPACEDIM == 3
          }
#endif

          coarFluxReg(iv, a_coarVar) *= invFinePerCoar;
        };

        auto irregularKernel = [&](const FaceIndex& face) -> void {
          coarFlux(face, a_coarVar) = 0.0;

          const FaceStencil& stencil = coarseningStencils(face, 0);
          for (int i = 0; i < stencil.size(); i++) {
            const FaceIndex& fineFace   = stencil.face(i);
            const Real       fineWeight = stencil.weight(i);

            coarFlux(face, a_coarVar) += fineWeight * fineFlux(fineFace, a_fineVar);
          }
        };

        // Kernel regions
        CH_START(t1);
        const Box     sideBox = adjCellBox(coarCellBox, dir, sit(), -1);
        const Box     faceBox = surroundingNodes(sideBox, dir);
        FaceIterator& faceIt  = m_fluxCoarseningRegions[din][std::make_pair(dir, sit())];
        CH_STOP(t1);

        CH_START(t2);
        BoxLoops::loop(faceBox, regularKernel);
        CH_STOP(t2);

        CH_START(t3);
        BoxLoops::loop(faceIt, irregularKernel);
        CH_STOP(t3);
      }
    }
  }
}

void
EBReflux::refluxIntoCoarse(LevelData<EBCellFAB>&       a_Lphi,
                           const LevelData<EBFluxFAB>& a_oldFluxes,
                           const LevelData<EBFluxFAB>& a_newFluxes,
                           const int                   a_phiVar,
                           const int                   a_oldFluxVar,
                           const int                   a_newFluxVar,
                           const Real                  a_scaleCoarFlux,
                           const Real                  a_scaleFineFlux) const noexcept
{
  CH_TIMERS("EBReflux::refluxIntoCoarse");
  CH_TIMER("EBReflux::refluxIntoCoarse::regular_cells", t1);
  CH_TIMER("EBReflux::refluxIntoCoarse::irregular_cells", t2);

  CH_assert(m_isDefined);
  CH_assert(a_Lphi.nComp() > a_phiVar);
  CH_assert(a_oldFluxes.nComp() > a_oldFluxVar);
  CH_assert(a_newFluxes.nComp() > a_newFluxVar);

  const DisjointBoxLayout& dbl   = m_eblg.getDBL();
  const DataIterator&      dit   = dbl.dataIterator();
  const EBISLayout&        ebisl = m_eblg.getEBISL();

  const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const Box      cellBox = dbl[din];
    const EBISBox& ebisBox = ebisl[din];

    EBCellFAB& Lphi    = a_Lphi[din];
    FArrayBox& LphiReg = Lphi.getFArrayBox();

    for (int dir = 0; dir < SpaceDim; dir++) {
      const EBFaceFAB& oldFlux = a_oldFluxes[din][dir];
      const EBFaceFAB& newFlux = a_newFluxes[din][dir];

      const FArrayBox& oldFluxReg = oldFlux.getFArrayBox();
      const FArrayBox& newFluxReg = newFlux.getFArrayBox();

      for (SideIterator sit; sit.ok(); ++sit) {
        const int     iHiLo = sign(flip(sit()));
        const IntVect shift = (sit() == Side::Lo) ? BASISV(dir) : IntVect::Zero;

        auto regularKernel = [&](const IntVect& iv) -> void {
          LphiReg(iv, a_phiVar) -= iHiLo * a_scaleCoarFlux * oldFluxReg(iv + shift, a_oldFluxVar);
          LphiReg(iv, a_phiVar) += iHiLo * a_scaleFineFlux * newFluxReg(iv + shift, a_newFluxVar);
        };

        auto irregularKernel = [&](const VolIndex& vof) -> void {
          const Vector<FaceIndex>& faces = ebisBox.getFaces(vof, dir, flip(sit()));

          for (int iface = 0; iface < faces.size(); iface++) {
            const FaceIndex& face     = faces[iface];
            const Real       faceArea = ebisBox.areaFrac(face);

            Lphi(vof, a_phiVar) -= iHiLo * faceArea * a_scaleCoarFlux * oldFlux(face, a_oldFluxVar);
            Lphi(vof, a_phiVar) += iHiLo * faceArea * a_scaleFineFlux * newFlux(face, a_newFluxVar);
          }
        };

        const DenseIntVectSet& regularCFIVS   = m_regularCoarseFineRegions[din].at(std::make_pair(dir, sit()));
        VoFIterator&           irregularCFIVS = m_irregularCoarseFineRegions[din].at(std::make_pair(dir, sit()));

        CH_START(t1);
        BoxLoops::loop(regularCFIVS, regularKernel);
        CH_STOP(t1);

        CH_START(t1);
        BoxLoops::loop(irregularCFIVS, irregularKernel);
        CH_STOP(t1);
      }
    }
  }
}

#include <CD_NamespaceFooter.H>
