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
#include <EBAlias.H>

// Our includes
#include <CD_EBCoarseFineParticleMesh.H>
#include <CD_EBAddOp.H>
#include <CD_DataOps.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

constexpr int EBCoarseFineParticleMesh::m_comp;
constexpr int EBCoarseFineParticleMesh::m_nComp;

EBCoarseFineParticleMesh::EBCoarseFineParticleMesh()
{
  CH_TIME("EBCoarseFineParticleMesh::EBCoarseFineParticleMesh");

  m_isDefined = false;
}

EBCoarseFineParticleMesh::EBCoarseFineParticleMesh(const EBLevelGrid& a_eblgCoar,
                                                   const EBLevelGrid& a_eblgFine,
                                                   const int          a_refRat,
                                                   const IntVect      a_ghost)
{
  CH_TIME("EBCoarseFineParticleMesh::EBCoarseFineParticleMesh");

  this->define(a_eblgCoar, a_eblgFine, a_refRat, a_ghost);
}

EBCoarseFineParticleMesh::~EBCoarseFineParticleMesh()
{
  CH_TIME("EBCoarseFineParticleMesh::~EBCoarseFineParticleMesh");
}

void
EBCoarseFineParticleMesh::define(const EBLevelGrid& a_eblgCoar,
                                 const EBLevelGrid& a_eblgFine,
                                 const int          a_refRat,
                                 const IntVect      a_ghost)
{
  CH_TIME("EBCoarseFineParticleMesh::define");

  CH_assert(a_refRat % 2 == 0);
  CH_assert(a_refRat >= 2);

  m_eblgCoar = a_eblgCoar;
  m_eblgFine = a_eblgFine;
  m_refRat   = a_refRat;
  m_ghost    = a_ghost;

  // Make the coarsened and refine grids.
  coarsen(m_eblgCoFi, m_eblgFine, m_refRat);
  refine(m_eblgFiCo, m_eblgCoar, m_refRat);

  // Create various buffers. The fine data should just be a raw copy.
  m_bufferCoFi.define(m_eblgCoFi.getDBL(), m_nComp, m_ghost, EBCellFactory(m_eblgCoFi.getEBISL()));
  m_bufferFiCo.define(m_eblgFiCo.getDBL(), m_nComp, m_ghost, EBCellFactory(m_eblgFiCo.getEBISL()));
  m_bufferFine.define(m_eblgFine.getDBL(), m_nComp, m_ghost, EBCellFactory(m_eblgFine.getEBISL()));

  // Define Copiers. Note that when we call ghostDefine we are actually adding from ghost+valid cells in the source to valid+ghost cells
  // in the destination.
  //
  // On the other hand, when we call define we are adding from the valid region in source to valid+ghost in the destination. These things
  // matter.

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

  // Define VoF iterators
  this->defineVoFIterators();

  m_isDefined = true;
}

void
EBCoarseFineParticleMesh::defineVoFIterators()
{
  CH_TIME("EBCoarseFineParticleMesh::defineVoFIterators");

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

  for (DataIterator dit(dblFine); dit.ok(); ++dit) {
    const Box&     cellBoxFine = dblFine[dit()];
    const EBISBox& ebisBoxFine = ebislFine[dit()];
    const EBGraph& ebgraphFine = ebisBoxFine.getEBGraph();

    const EBISBox& ebisBoxCoFi = ebislCoFi[dit()];
    const EBGraph& ebgraphCoFi = ebisBoxCoFi.getEBGraph();

    Box grownBoxFine = grow(cellBoxFine, m_ghost);
    grownBoxFine &= domainFine;

    // On the fine grid I only want the ghost cells that are irregular.
    IntVectSet irregIVSFine = ebisBoxFine.getIrregIVS(grownBoxFine);
    irregIVSFine -= cellBoxFine;
    irregIVSFine &= domainFine;
    m_vofIterFineGhosts[dit()].define(irregIVSFine, ebgraphFine);

    // On the coarse grid I want coarsenings of the irregular ghost cells on the fine level.
    IntVectSet ivsCoFi = coarsen(irregIVSFine, m_refRat);
    ivsCoFi &= domainCoar;
    m_vofIterCoFiGhosts[dit()].define(ivsCoFi, ebgraphCoFi);
  }

  for (DataIterator dit(dblCoar); dit.ok(); ++dit) {
    const Box&       cellBoxCoar  = dblCoar[dit()];
    const EBISBox&   ebisBoxCoar  = ebislCoar[dit()];
    const EBGraph&   ebGraphCoar  = ebisBoxCoar.getEBGraph();
    const IntVectSet irregIVSCoar = ebisBoxCoar.getIrregIVS(cellBoxCoar);

    m_vofIterCoar[dit()].define(irregIVSCoar, ebGraphCoar);
  }
}

void
EBCoarseFineParticleMesh::addFineGhostsToCoarse(LevelData<EBCellFAB>&       a_coarData,
                                                const LevelData<EBCellFAB>& a_fineData) const
{
  CH_TIME("EBCoarseFineParticleMesh::addFineGhostsToCoarse");

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() == 1);
  CH_assert(a_fineData.nComp() == 1);
  CH_assert(m_refRat % 2 == 0);
  CH_assert(a_coarData.ghostVect() == m_ghost);
  CH_assert(a_fineData.ghostVect() == m_ghost);

  // TLDR: This routine will take the ghost cells that are in a_fineData and add them to the coarse level. We use buffers for this,
  //       copying the data from a_fineData to our nifty m_bufferFine buffer holder. The ghost cells in that scratch data are
  //       coarsened onto yet another buffer (m_bufferCoFi). Finally, we add the contents in that buffer to the coarse data.

  const DisjointBoxLayout& dblFine    = m_eblgFine.getDBL();
  const ProblemDomain&     domainFine = m_eblgFine.getDomain();
  const EBISLayout&        ebislFine  = m_eblgFine.getEBISL();

  const DisjointBoxLayout& dblCoFi   = m_eblgCoFi.getDBL();
  const EBISLayout&        ebislCoFi = m_eblgCoFi.getEBISL();

  // Copy the fine data to scratch and reset the interior cells. We do this by copying everything to the fine scratch data,
  // and then set all the valid data in each box to zero. The exchange() operation will then take care of ghost cells that.
  // overlap with valid regions in different boxes (we could use a NeighborIterator to achieve the same). Note that this is critical
  // for the result because we essentially end up setting all valid data to zero so we don't accidentally add mass to invalid region
  // of the coarse grid (i.e., the region underneath the fine grid).
  a_fineData.localCopyTo(m_bufferFine);
  for (DataIterator dit(dblFine); dit.ok(); ++dit) {
    const Box& box = dblFine[dit()];
    m_bufferFine[dit()].setVal(0.0, box, 0, m_nComp);
  }
  m_bufferFine.exchange();

  // Coarsen the fine grid data. We do this by AVERAGING the fine-grid data onto the coarse grid. This actually involves entire patches
  // and not just individual ghost cells regions, but that's ok because the rest of the data (in the valid fine regions) are set to zero above,
  // so the coarse data underneath those regions will be zero, also.
  const Real factor = 1. / pow(m_refRat, SpaceDim);
  for (DataIterator dit(dblFine); dit.ok(); ++dit) {
    EBCellFAB&       coFiData = m_bufferCoFi[dit()];
    const EBCellFAB& fineData = m_bufferFine[dit()];

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
    VoFIterator& vofitCoar = m_vofIterCoFiGhosts[dit()];
    VoFIterator& vofitFine = m_vofIterFineGhosts[dit()];

    // Set the value in the coarse cell.
    auto setCoarData = [&](const VolIndex& fineVof) -> void {
      const VolIndex coarVof = ebislFine.coarsen(fineVof, m_refRat, dit());

      for (int comp = 0; comp < m_nComp; comp++) {
        coFiData(coarVof, comp) = 0.0;
      }
    };

    // Add mass from the fine vof to the coarse vof.
    auto addFineToCoar = [&](const VolIndex& fineVof) -> void {
      const VolIndex& coarVof = ebislFine.coarsen(fineVof, m_refRat, dit());
      for (int comp = 0; comp < m_nComp; comp++) {
        coFiData(coarVof, comp) += fineData(fineVof, comp);
      }
    };

    // Scale the value in the coarse cell.
    auto scaleCoar = [&](const VolIndex& coarVof) -> void {
      const Vector<VolIndex> fineVofs = ebislCoFi.refine(coarVof, m_refRat, dit());

      for (int comp = 0; comp < m_nComp; comp++) {
        coFiData(coarVof, comp) *= 1. / fineVofs.size();
      }
    };

    // Execute kernels.
    BoxLoops::loop(vofitFine, setCoarData);
    BoxLoops::loop(vofitFine, addFineToCoar);
    BoxLoops::loop(vofitCoar, scaleCoar);
  }

  // At this point we have coarsened the data in the ghost regions around the fine grid onto m_bufferCoFi. We should be table to take that data and add
  // it directly to
  const Interval interv(0, m_nComp - 1);

  m_bufferCoFi.copyTo(interv, a_coarData, interv, m_copierCoFiToCoarIncludeGhosts, EBAddOp());
}

void
EBCoarseFineParticleMesh::addFiCoDataToFine(LevelData<EBCellFAB>&       a_fineData,
                                            const LevelData<EBCellFAB>& a_fiCoData) const
{
  CH_TIME("EBCoarseFineParticleMesh::addFiCoDataToFine");

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
                                                 const LevelData<EBCellFAB>& a_coarData) const
{
  CH_TIME("EBCoarseFineParticleMesh::addInvalidCoarseToFine");

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

  for (DataIterator dit(dblCoar); dit.ok(); ++dit) {
    EBCellFAB&       fiCoData = m_bufferFiCo[dit()];
    const EBCellFAB& coarData = a_coarData[dit()];
    const Box        coarBox  = dblCoar[dit()];

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
    VoFIterator& vofit = m_vofIterCoar[dit()];

    auto kernel = [&](const VolIndex& coarVoF) -> void {
      const Vector<VolIndex>& fineVoFs = ebislCoar.refine(coarVoF, m_refRat, dit());

      for (int ivof = 0; ivof < fineVoFs.size(); ivof++) {
        const VolIndex& fineVoF   = fineVoFs[ivof];
        fiCoData(fineVoF, m_comp) = coarData(coarVoF, m_comp);
      }
    };

    BoxLoops::loop(vofit, kernel);
  }

  // Finally, add the data to the valid region on the fine grid.
  const Interval interv(0, m_nComp - 1);
  m_bufferFiCo.copyTo(interv, a_fineData, interv, m_copierFiCoToFineNoGhosts, EBAddOp());
}

LevelData<EBCellFAB>&
EBCoarseFineParticleMesh::getFiCoBuffer() const
{
  CH_TIME("EBCoarseFineParticleMesh::getFiCoBuffer");

  return m_bufferFiCo;
}

const EBLevelGrid&
EBCoarseFineParticleMesh::getEblgFiCo() const
{
  CH_TIME("EBCoarseFineParticleMesh::getEblgFiCo");

  return m_eblgFiCo;
}

#include <CD_NamespaceFooter.H>
