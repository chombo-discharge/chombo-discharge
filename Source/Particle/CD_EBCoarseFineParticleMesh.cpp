/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBCoarseFineParticleMesh.cpp
  @brief  Implementation of CD_EBCoarseFineParticleMesh.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <EBCellFactory.H>
#include <EBAverageF_F.H>
#include <ParmParse.H>
#include <EBAlias.H>

// Our includes
#include <CD_EBCoarseFineParticleMesh.H>
#include <CD_EBAddOp.H>
#include <CD_DataOps.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

constexpr int EBCoarseFineParticleMesh::m_comp;
constexpr int EBCoarseFineParticleMesh::m_nComp;

EBCoarseFineParticleMesh::EBCoarseFineParticleMesh() noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::EBCoarseFineParticleMesh");

  m_isDefined = false;
  m_verbose   = false;
}

EBCoarseFineParticleMesh::EBCoarseFineParticleMesh(const EBLevelGrid& a_eblgCoar,
                                                   const EBLevelGrid& a_eblgFine,
                                                   const int          a_refRat,
                                                   const IntVect      a_ghost) noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::EBCoarseFineParticleMesh");

  m_verbose = false;

  this->define(a_eblgCoar, a_eblgFine, a_refRat, a_ghost);
}

EBCoarseFineParticleMesh::~EBCoarseFineParticleMesh() noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::~EBCoarseFineParticleMesh");
}

void
EBCoarseFineParticleMesh::define(const EBLevelGrid& a_eblgCoar,
                                 const EBLevelGrid& a_eblgFine,
                                 const int          a_refRat,
                                 const IntVect      a_ghost) noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::define");

  ParmParse pp("EBCoarseFineParticleMesh");
  pp.query("verbose", m_verbose);

  if (m_verbose) {
    pout() << "EBCoarseFineParticleMesh::define" << endl;
  }

  CH_assert(a_refRat % 2 == 0);
  CH_assert(a_refRat >= 2);

  m_eblgCoar = a_eblgCoar;
  m_eblgFine = a_eblgFine;
  m_refRat   = a_refRat;
  m_ghost    = a_ghost;

  // Make the coarsened and refine grids.
  coarsen(m_eblgCoFi, m_eblgFine, m_refRat);
  refine(m_eblgFiCo, m_eblgCoar, m_refRat);

  // valid+ghost -> valid
  m_copierCoFiToCoarIncludeGhosts.ghostDefine(m_eblgCoFi.getDBL(),
                                              m_eblgCoar.getDBL(),
                                              m_eblgCoar.getDomain(),
                                              m_ghost);
  // valid+ghost -> valid
  m_copierFiCoToFineIncludeGhosts.ghostDefine(m_eblgFiCo.getDBL(),
                                              m_eblgFine.getDBL(),
                                              m_eblgFine.getDomain(),
                                              m_ghost);

  // valid -> valid+ghost
  m_copierFiCoToFineNoGhosts.define(m_eblgFiCo.getDBL(), m_eblgFine.getDBL(), m_eblgFine.getDomain(), m_ghost);

  // valid+ghost -> valid on the m_eblgFICo only.
  m_copierFiCoToFiCo.define(m_eblgFiCo.getDBL(), m_eblgFiCo.getDBL(), m_eblgFine.getDomain(), m_ghost, true);
  m_copierFiCoToFiCo.reverse();

  // Define VoF iterators and stencils.
  this->defineVoFIterators();
  this->defineStencils();
  this->defineBufferFiCo();

  m_isDefined = true;
}

void
EBCoarseFineParticleMesh::defineVoFIterators() noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::defineVoFIterators");
  if (m_verbose) {
    pout() << "EBCoarseFineParticleMesh::defineVoFIterators" << endl;
  }

  const DisjointBoxLayout& dblCoar    = m_eblgCoar.getDBL();
  const ProblemDomain&     domainCoar = m_eblgCoar.getDomain();
  const EBISLayout&        ebislCoar  = m_eblgCoar.getEBISL();

  const DisjointBoxLayout& dblFine    = m_eblgFine.getDBL();
  const ProblemDomain&     domainFine = m_eblgFine.getDomain();
  const EBISLayout&        ebislFine  = m_eblgFine.getEBISL();

  const DisjointBoxLayout& dblCoFi   = m_eblgCoFi.getDBL();
  const EBISLayout&        ebislCoFi = m_eblgCoFi.getEBISL();

  m_vofIterFineGhosts.define(dblFine);
  m_vofIterCoFiGhosts.define(dblCoFi);
  m_vofIterCoar.define(dblCoar);

  const DataIterator& ditFine = dblFine.dataIterator();
  const DataIterator& ditCoar = dblCoar.dataIterator();

  const int nboxFine = ditFine.size();
  const int nboxCoar = ditCoar.size();

#pragma omp parallel
  {
#pragma omp for schedule(runtime)
    for (int mybox = 0; mybox < nboxFine; mybox++) {
      const DataIndex& din = ditFine[mybox];

      const Box&     cellBoxFine = dblFine[din];
      const EBISBox& ebisBoxFine = ebislFine[din];
      const EBGraph& ebgraphFine = ebisBoxFine.getEBGraph();

      const EBISBox& ebisBoxCoFi = ebislCoFi[din];
      const EBGraph& ebgraphCoFi = ebisBoxCoFi.getEBGraph();

      Box grownBoxFine = grow(cellBoxFine, m_ghost);
      grownBoxFine &= domainFine;

      // On the fine grid I only want the ghost cells that are irregular.
      IntVectSet irregIVSFine = ebisBoxFine.getIrregIVS(grownBoxFine);
      irregIVSFine -= cellBoxFine;
      irregIVSFine &= domainFine;
      m_vofIterFineGhosts[din].define(irregIVSFine, ebgraphFine);

      // On the coarse grid I want coarsenings of the irregular ghost cells on the fine level.
      IntVectSet ivsCoFi = coarsen(irregIVSFine, m_refRat);
      ivsCoFi &= domainCoar;
      m_vofIterCoFiGhosts[din].define(ivsCoFi, ebgraphCoFi);
    }

#pragma omp for schedule(runtime)
    for (int mybox = 0; mybox < nboxCoar; mybox++) {
      const DataIndex& din = ditCoar[mybox];

      const Box&       cellBoxCoar  = dblCoar[din];
      const EBISBox&   ebisBoxCoar  = ebislCoar[din];
      const EBGraph&   ebGraphCoar  = ebisBoxCoar.getEBGraph();
      const IntVectSet irregIVSCoar = ebisBoxCoar.getIrregIVS(cellBoxCoar);

      m_vofIterCoar[din].define(irregIVSCoar, ebGraphCoar);
    }
  }
}

void
EBCoarseFineParticleMesh::defineStencils() noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::defineStencils");
  if(m_verbose) {
    pout() << "EBCoarseFineParticleMesh::defineStencils" << endl;
  }

  const DisjointBoxLayout& dblCoar = m_eblgCoar.getDBL();
  const DisjointBoxLayout& dblFiCo = m_eblgFiCo.getDBL();

  const DataIterator& ditCoar = dblCoar.dataIterator();

  const EBISLayout& ebislCoar = m_eblgCoar.getEBISL();
  const EBISLayout& ebislFiCo = m_eblgFiCo.getEBISL();

  const Real dxCoar   = 1.0;
  const Real dxFine   = dxCoar / m_refRat;
  const Real dxFactor = std::pow(dxFine / dxCoar, SpaceDim);
  const int  nbox     = ditCoar.size();

  m_cellArithmeticStencils.define(dblCoar);
  m_cellConservativeStencils.define(dblCoar);

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex&  din         = ditCoar[mybox];
    const EBISBox&    ebisBoxCoar = ebislCoar[din];
    const EBISBox&    ebisBoxFiCo = ebislFiCo[din];
    const EBGraph&    ebGraphCoar = ebisBoxCoar.getEBGraph();
    const Box&        cellBox     = dblCoar[din];
    const IntVectSet& irregCoar   = ebisBoxCoar.getIrregIVS(cellBox);

    BaseIVFAB<VoFStencil>& arithmeticStencils   = m_cellArithmeticStencils[din];
    BaseIVFAB<VoFStencil>& conservativeStencils = m_cellConservativeStencils[din];

    arithmeticStencils.define(irregCoar, ebGraphCoar, 1);
    conservativeStencils.define(irregCoar, ebGraphCoar, 1);

    auto buildStencils = [&](const VolIndex& coarVoF) -> void {
      const Real             kappaC      = ebisBoxCoar.volFrac(coarVoF);
      const Vector<VolIndex> fineVoFs    = ebislCoar.refine(coarVoF, m_refRat, din);
      const int              numFineVoFs = fineVoFs.size();

      VoFStencil& arithSten = arithmeticStencils(coarVoF, 0);
      VoFStencil& consSten  = conservativeStencils(coarVoF, 0);

      arithSten.clear();
      consSten.clear();

      if (numFineVoFs > 0) {

        for (int ifine = 0; ifine < numFineVoFs; ifine++) {
          const VolIndex& fineVoF = fineVoFs[ifine];
          const Real      kappaF  = ebisBoxFiCo.volFrac(fineVoF);

          arithSten.add(fineVoF, 1.0 / numFineVoFs);

          if (kappaC > 0.0) {
            consSten.add(fineVoF, kappaF * dxFactor / kappaC);
          }
          else {
            consSten.add(fineVoF, 1.0 / numFineVoFs);
          }
        }
      }
    };

    BoxLoops::loop(m_vofIterCoar[din], buildStencils);
  }
}

void
EBCoarseFineParticleMesh::defineBufferFiCo() noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::defineBufferFiCo");
  if (m_verbose) {
    pout() << "EBCoarseFineParticleMesh::defineBufferFiCo" << endl;
  }

  const DisjointBoxLayout& dblFiCo   = m_eblgFiCo.getDBL();
  const EBISLayout&        ebislFiCo = m_eblgFiCo.getEBISL();

  m_bufferFiCoReal.define(dblFiCo, 1, m_ghost * IntVect::Unit, EBCellFactory(ebislFiCo));
  m_bufferFiCoRealVect.define(dblFiCo, SpaceDim, m_ghost * IntVect::Unit, EBCellFactory(ebislFiCo));
}

void
EBCoarseFineParticleMesh::addFineGhostsToCoarse(LevelData<EBCellFAB>&       a_coarData,
                                                const LevelData<EBCellFAB>& a_fineData) const noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::addFineGhostsToCoarse");
  if (m_verbose) {
    pout() << "EBCoarseFineParticleMesh::addFineGhostsToCoarse" << endl;
  }

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() == 1);
  CH_assert(a_fineData.nComp() == 1);
  CH_assert(m_refRat % 2 == 0);
  CH_assert(a_coarData.ghostVect() == m_ghost);
  CH_assert(a_fineData.ghostVect() == m_ghost);

  // TLDR: This routine will take the ghost cells that are in a_fineData and add them to the coarse level. We use buffers for this,
  //       copying the data from a_fineData to our nifty m_bufferFine buffer holder. The ghost cells in that scratch data are
  //       coarsened onto yet another buffer. Finally, we add the contents in that buffer to the coarse data.

  const DisjointBoxLayout& dblFine    = m_eblgFine.getDBL();
  const ProblemDomain&     domainFine = m_eblgFine.getDomain();
  const EBISLayout&        ebislFine  = m_eblgFine.getEBISL();

  const DisjointBoxLayout& dblCoFi   = m_eblgCoFi.getDBL();
  const EBISLayout&        ebislCoFi = m_eblgCoFi.getEBISL();

  const DataIterator& ditFine = dblFine.dataIterator();
  const DataIterator& ditCoFi = dblCoFi.dataIterator();

  const int nboxFine = ditFine.size();
  const int nboxCoFi = ditCoFi.size();

  const Real factor = 1. / pow(m_refRat, SpaceDim);

  LevelData<EBCellFAB> bufferFine(m_eblgFine.getDBL(), m_nComp, m_ghost, EBCellFactory(m_eblgFine.getEBISL()));
  LevelData<EBCellFAB> bufferCoFi(m_eblgCoFi.getDBL(), m_nComp, m_ghost, EBCellFactory(m_eblgCoFi.getEBISL()));

  // Copy the fine data to scratch and reset the interior cells. We do this by copying everything to the fine scratch data,
  // and then set all the valid data in each box to zero. The exchange() operation will then take care of ghost cells that.
  // overlap with valid regions in different boxes (we could use a NeighborIterator to achieve the same). Note that this is critical
  // for the result because we essentially end up setting all valid data to zero so we don't accidentally add mass to invalid region
  // of the coarse grid (i.e., the region underneath the fine grid).
  a_fineData.localCopyTo(bufferFine);
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nboxFine; mybox++) {
    const DataIndex& din = ditFine[mybox];

    const Box& box = dblFine[din];

    bufferFine[din].setVal(0.0, box, 0, m_nComp);
  }
  bufferFine.exchange();

  // Coarsen the fine grid data. We do this by AVERAGING the fine-grid data onto the coarse grid. This actually involves entire patches
  // and not just individual ghost cells regions, but that's ok because the rest of the data (in the valid fine regions) are set to zero above,
  // so the coarse data underneath those regions will be zero, also.
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nboxFine; mybox++) {
    const DataIndex& din = ditFine[mybox];

    EBCellFAB&       coFiData = bufferCoFi[din];
    const EBCellFAB& fineData = bufferFine[din];

    coFiData.setVal(0.0);

    // Do all the regular cells.
    FArrayBox&       coFiDataReg = coFiData.getFArrayBox();
    const FArrayBox& fineDataReg = fineData.getFArrayBox();

    const Box fineBox = fineDataReg.box() & domainFine;
    for (int comp = 0; comp < m_nComp; comp++) {

      auto regularKernel = [&](const IntVect& ivFine) -> void {
        const IntVect ivCoar = coarsen(ivFine, m_refRat);

        coFiDataReg(ivCoar, comp) += fineDataReg(ivFine, comp) * factor;
      };

      BoxLoops::loop(fineBox, regularKernel);
    }

    // Now do the irregular cells. First reset the coarse cell values and then increment with the fine-grid values.
    VoFIterator& vofitCoar = m_vofIterCoFiGhosts[din];
    VoFIterator& vofitFine = m_vofIterFineGhosts[din];

    // Set the value in the coarse cell.
    auto setCoarData = [&](const VolIndex& fineVof) -> void {
      const VolIndex coarVof = ebislFine.coarsen(fineVof, m_refRat, din);

      for (int comp = 0; comp < m_nComp; comp++) {
        coFiData(coarVof, comp) = 0.0;
      }
    };

    // Add mass from the fine vof to the coarse vof.
    auto addFineToCoar = [&](const VolIndex& fineVof) -> void {
      const VolIndex& coarVof = ebislFine.coarsen(fineVof, m_refRat, din);
      for (int comp = 0; comp < m_nComp; comp++) {
        coFiData(coarVof, comp) += fineData(fineVof, comp);
      }
    };

    // Scale the value in the coarse cell.
    auto scaleCoar = [&](const VolIndex& coarVof) -> void {
      const Vector<VolIndex> fineVofs = ebislCoFi.refine(coarVof, m_refRat, din);

      for (int comp = 0; comp < m_nComp; comp++) {
        coFiData(coarVof, comp) *= 1. / fineVofs.size();
      }
    };

    // Execute kernels.
    BoxLoops::loop(vofitFine, setCoarData);
    BoxLoops::loop(vofitFine, addFineToCoar);
    BoxLoops::loop(vofitCoar, scaleCoar);
  }

  const Interval interv(0, m_nComp - 1);

  bufferCoFi.copyTo(interv, a_coarData, interv, m_copierCoFiToCoarIncludeGhosts, EBAddOp());
}

const EBLevelGrid&
EBCoarseFineParticleMesh::getEblgFiCo() const
{
  CH_TIME("EBCoarseFineParticleMesh::getEblgFiCo");
  if (m_verbose) {
    pout() << "EBCoarseFineParticleMesh::getEblgFiCo" << endl;
  }

  return m_eblgFiCo;
}

template <>
LevelData<EBCellFAB>&
EBCoarseFineParticleMesh::getBufferFiCo<1>() const noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::getBufferFiCo<1>");
  if (m_verbose) {
    pout() << "EBCoarseFineParticleMesh::getBufferFiCo<1>" << endl;
  }

  return m_bufferFiCoReal;
}

template <>
LevelData<EBCellFAB>&
EBCoarseFineParticleMesh::getBufferFiCo<SpaceDim>() const noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::getBufferFiCo<SpaceDim>");
  if (m_verbose) {
    pout() << "EBCoarseFineParticleMesh::getBufferFiCo<SpaceDim>" << endl;
  }

  return m_bufferFiCoRealVect;
}

void
EBCoarseFineParticleMesh::addFiCoDataToFine(LevelData<EBCellFAB>&       a_fineData,
                                            const LevelData<EBCellFAB>& a_fiCoData) const noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::addFiCoDataToFine");
  if (m_verbose) {
    pout() << "EBCoarseFineParticleMesh::addFiCoDataToFine" << endl;
  }

  CH_assert(m_isDefined);
  CH_assert(a_fineData.nComp() == 1);
  CH_assert(a_fiCoData.nComp() == 1);
  CH_assert(a_fineData.ghostVect() == m_ghost);
  CH_assert(a_fiCoData.ghostVect() == m_ghost);

  const Interval interv(0, m_nComp - 1);

  a_fiCoData.copyTo(interv, a_fineData, interv, m_copierFiCoToFineIncludeGhosts, EBAddOp());
}

void
EBCoarseFineParticleMesh::addInvalidCoarseToFine(LevelData<EBCellFAB>&       a_fineData,
                                                 const LevelData<EBCellFAB>& a_coarData) const noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::addInvalidCoarseToFine");
  if (m_verbose) {
    pout() << "EBCoarseFineParticleMesh::addInvalidCoarseToFine" << endl;
  }

  CH_assert(m_isDefined);
  CH_assert(a_fineData.nComp() == 1);
  CH_assert(a_coarData.nComp() == 1);
  CH_assert(a_fineData.ghostVect() == m_ghost);
  CH_assert(a_coarData.ghostVect() == m_ghost);

  // TLDR: This routine performs a piecewise constant interpolation of the coarse data to the fine grid. We do this by going through the coarse-grid data
  //       and piecewise interpolating the result to the fine grid (using a buffer). After that, we add the contents in the buffer to the fine level.
  //
  //       The data-motion plan for this is to add the valid+ghost cells in the interpolated fine-grid data to the valid region on the fine grid. Note that
  //       the function signature indicates that we only add invalid coarse data, we do run kernels over all data. The addition is done only at the end in
  //       copyTo.
  const DisjointBoxLayout& dblCoar   = m_eblgCoar.getDBL();
  const EBISLayout&        ebislCoar = m_eblgCoar.getEBISL();

  const DisjointBoxLayout& dblFiCo   = m_eblgFiCo.getDBL();
  const EBISLayout&        ebislFiCo = m_eblgFiCo.getEBISL();

  LevelData<EBCellFAB> bufferFiCo(m_eblgFiCo.getDBL(), m_nComp, m_ghost, EBCellFactory(m_eblgFiCo.getEBISL()));

  const DataIterator dit = dblCoar.dataIterator();

  const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    EBCellFAB&       fiCoData = bufferFiCo[din];
    const EBCellFAB& coarData = a_coarData[din];
    const Box        coarBox  = dblCoar[din];

    fiCoData.setVal(0.0);

    // Regular cells.
    BaseFab<Real>&       fiCoDataReg = fiCoData.getSingleValuedFAB();
    const BaseFab<Real>& coarDataReg = coarData.getSingleValuedFAB();

    // This is the regular kernel -- it sets the fine data equal to the coarse data.
    auto regularKernel = [&](const IntVect& ivCoar) -> void {
      // May seem weird, but we are running nested box loops here. The second kernel runs over the refined box Box(IntVect::Zero, m_refRat*IntVect::Unit),
      // which gives the number of fine-grid cells that lie on top of a coarse grid cell. So, we are iterating over that box and figuring out which cells in
      // that box correspond to which fine cells (we just need to increment by m_refRat * ivCoar).
      auto fineKernel = [&](const IntVect& iv) {
        const IntVect ivFine = m_refRat * ivCoar + iv;

        fiCoDataReg(ivFine, m_comp) = coarDataReg(ivCoar, m_comp);
      };

      BoxLoops::loop(Box(IntVect::Zero, m_refRat * IntVect::Unit), fineKernel);
    };

    // Execute kernel over the entire coarse-grid patch.
    BoxLoops::loop(coarBox, regularKernel);

    // Now do the irregular cells. Here, we loop over all the coarse cells (including ghosts) and set the value in the
    // fine cells to be the same as the value in the underlying coarse cell.
    VoFIterator& vofit = m_vofIterCoar[din];

    auto kernel = [&](const VolIndex& coarVoF) -> void {
      const Vector<VolIndex>& fineVoFs = ebislCoar.refine(coarVoF, m_refRat, din);

      for (int ivof = 0; ivof < fineVoFs.size(); ivof++) {
        const VolIndex& fineVoF   = fineVoFs[ivof];
        fiCoData(fineVoF, m_comp) = coarData(coarVoF, m_comp);
      }
    };

    BoxLoops::loop(vofit, kernel);
  }

  // Finally, add the data to the valid region on the fine grid.
  const Interval interv(0, m_nComp - 1);
  bufferFiCo.copyTo(interv, a_fineData, interv, m_copierFiCoToFineNoGhosts, EBAddOp());
}

void
EBCoarseFineParticleMesh::exchangeAndAddFiCoData(LevelData<EBCellFAB>& a_fiCoData) const noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::addInvalidCoarseToFine");
  if (m_verbose) {
    pout() << "EBCoarseFineParticleMesh::addInvalidCoarseToFine" << endl;
  }

  CH_assert(m_isDefined);
  CH_assert(a_fiCoData.ghostVect() == m_ghost);
  CH_assert(a_fiCoData.nComp() == 1);

  const Interval interv(0, m_nComp - 1);

  a_fiCoData.exchange(interv, m_copierFiCoToFiCo, EBAddOp());
}

void
EBCoarseFineParticleMesh::restrictAndAddFiCoDataToCoar(LevelData<EBCellFAB>&                   a_coarData,
                                                       const LevelData<EBCellFAB>&             a_fiCoData,
                                                       const EBCoarseFineParticleMesh::Average a_average) const noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::restrictAndAddFiCoDataToCoar");
  if (m_verbose) {
    pout() << "EBCoarseFineParticleMesh::restrictAndAddFiCoDataToCoar" << endl;
  }

  const DisjointBoxLayout& dblCoar = m_eblgCoar.getDBL();
  const DisjointBoxLayout& dblFiCo = m_eblgFiCo.getDBL();

  CH_assert(m_isDefined);
  CH_assert(a_coarData.ghostVect() == m_ghost);
  CH_assert(a_fiCoData.ghostVect() == m_ghost);
  CH_assert(a_fiCoData.nComp() == 1);
  CH_assert(a_coarData.nComp() == 1);
  CH_assert(a_coarData.disjointBoxLayout() == dblCoar);
  CH_assert(a_fiCoData.disjointBoxLayout() == dblFiCo);

  const DataIterator& dit      = dblCoar.dataIterator();
  const int           numBoxes = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < numBoxes; mybox++) {
    const DataIndex& din = dit[mybox];

    EBCellFAB&       coarData = a_coarData[din];
    const EBCellFAB& fineData = a_fiCoData[din];

    switch (a_average) {
    case EBCoarseFineParticleMesh::Average::Arithmetic: {
      this->arithmeticAverageAndAdd(coarData, fineData, din);

      break;
    }
    case EBCoarseFineParticleMesh::Average::Conservative: {
      this->conservativeAverageAndAdd(coarData, fineData, din);

      break;
    }
    default: {
      MayDay::Abort("EBCoarseFineParticleMesh::restrictAndAddFiCoDataToCoar - logic bust in average method");

      break;
    }
    }
  }
}

void
EBCoarseFineParticleMesh::conservativeAverageAndAdd(EBCellFAB&       a_coarData,
                                                    const EBCellFAB& a_fineData,
                                                    const DataIndex& a_din) const noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::conservativeAverageAndAdd");
  if (m_verbose) {
    pout() << "EBCoarseFineParticleMesh::conservativeAverageAndAdd" << endl;
  }

  const Real dxCoar   = 1.0;
  const Real dxFine   = 1.0 / m_refRat;
  const Real dxFactor = std::pow(dxFine / dxCoar, SpaceDim);
  const Box  refiBox  = Box(IntVect::Zero, (m_refRat - 1) * IntVect::Unit);

  FArrayBox&       coarDataReg = a_coarData.getFArrayBox();
  const FArrayBox& fineDataReg = a_fineData.getFArrayBox();

  const BaseIVFAB<VoFStencil>& stencils = m_cellConservativeStencils[a_din];

  // Kernel for regular grid cells
  auto regularKernel = [&](const IntVect& iv) -> void {
    for (BoxIterator bit(refiBox); bit.ok(); ++bit) {
      coarDataReg(iv, 0) += fineDataReg(m_refRat * iv + bit(), 0) * dxFactor;
    }
  };

  // Kernel for irregular cells
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    const VoFStencil& stencil = stencils(vof, 0);

    for (int i = 0; i < stencil.size(); i++) {
      const VolIndex& ivof    = stencil.vof(i);
      const Real&     iweight = stencil.weight(i);

      a_coarData(vof, 0) += iweight * a_fineData(ivof, 0);
    }
  };

  BoxLoops::loop(m_eblgCoar.getDBL()[a_din], regularKernel);
  BoxLoops::loop(m_vofIterCoar[a_din], irregularKernel);
}

void
EBCoarseFineParticleMesh::arithmeticAverageAndAdd(EBCellFAB&       a_coarData,
                                                  const EBCellFAB& a_fineData,
                                                  const DataIndex& a_din) const noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::arithmeticAverageAndAdd");
  if (m_verbose) {
    pout() << "EBCoarseFineParticleMesh::arithmeticAverageAndAdd" << endl;
  }

  const Real dxCoar   = 1.0;
  const Real dxFine   = 1.0 / m_refRat;
  const Real dxFactor = std::pow(dxFine / dxCoar, SpaceDim);
  const Box  refiBox  = Box(IntVect::Zero, (m_refRat - 1) * IntVect::Unit);

  FArrayBox&       coarDataReg = a_coarData.getFArrayBox();
  const FArrayBox& fineDataReg = a_fineData.getFArrayBox();

  const BaseIVFAB<VoFStencil>& stencils = m_cellArithmeticStencils[a_din];

  // Kernel for regular grid cells
  auto regularKernel = [&](const IntVect& iv) -> void {
    for (BoxIterator bit(refiBox); bit.ok(); ++bit) {
      coarDataReg(iv, 0) += fineDataReg(m_refRat * iv + bit(), 0) * dxFactor;
    }
  };

  // Kernel for irregular cells
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    const VoFStencil& stencil = stencils(vof, 0);

    for (int i = 0; i < stencil.size(); i++) {
      const VolIndex& ivof    = stencil.vof(i);
      const Real&     iweight = stencil.weight(i);

      a_coarData(vof, 0) += iweight * a_fineData(ivof, 0);
    }
  };

  BoxLoops::loop(m_eblgCoar.getDBL()[a_din], regularKernel);
  BoxLoops::loop(m_vofIterCoar[a_din], irregularKernel);
}

#include <CD_NamespaceFooter.H>
