/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_Realm.cpp
  @brief  Implementation of CD_Realm.H
  @author Robert Marskar
*/

// Std includes
#include <cmath>

// Chombo includes
#include <ParmParse.H>
#include <NeighborIterator.H>
#include <BoxIterator.H>
#include <IntVectSet.H>
#include <DenseIntVectSet.H>
#include <Copier.H>

// Our includes
#include <CD_Realm.H>
#include <CD_BoxLoops.H>
#include <CD_OpenMP.H>
#include <CD_MemoryReport.H>
#include <CD_NamespaceHeader.H>

const std::string Realm::Primal = "primal";
const std::string Realm::primal = "primal";

Realm::Realm() : m_isDefined(false), m_verbosity(-1)
{

  ParmParse pp("Realm");

  pp.query("verbosity", m_verbosity);

  // Just empty points until define() is called
  m_realms.emplace(phase::gas, RefCountedPtr<PhaseRealm>(new PhaseRealm()));
  m_realms.emplace(phase::solid, RefCountedPtr<PhaseRealm>(new PhaseRealm()));
}

Realm::~Realm() = default;

void
Realm::define(const Vector<DisjointBoxLayout>&                           a_grids,
              const Vector<ProblemDomain>&                               a_domains,
              const Vector<int>&                                         a_refRat,
              const Vector<Real>&                                        a_dx,
              const RealVect&                                            a_probLo,
              const int                                                  a_finestLevel,
              const int                                                  a_minBlockSize,
              const int                                                  a_ebGhost,
              const int                                                  a_numGhost,
              const int                                                  a_lsfGhost,
              const int                                                  a_redistRad,
              const int                                                  a_mgInterpOrder,
              const int                                                  a_mgInterpRadius,
              const int                                                  a_mgInterpWeight,
              const CellCentroidInterpolation::Type                      a_centroidInterpType,
              const EBCentroidInterpolation::Type                        a_ebInterpType,
              const std::map<phase::which_phase, RefCountedPtr<BaseIF>>& a_baseif,
              const RefCountedPtr<MultiFluidIndexSpace>&                 a_mfis)
{
  CH_TIME("Realm::define");
  if (m_verbosity > 5) {
    pout() << "Realm::define" << endl;
  }

  m_grids                = a_grids;
  m_domains              = a_domains;
  m_refinementRatios     = a_refRat;
  m_dx                   = a_dx;
  m_probLo               = a_probLo;
  m_finestLevel          = a_finestLevel;
  m_minBlockSize         = a_minBlockSize;
  m_numGhost             = a_numGhost;
  m_baseif               = a_baseif;
  m_multifluidIndexSpace = a_mfis;

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  m_realms[phase::gas]->define(a_grids,
                               a_domains,
                               a_refRat,
                               a_dx,
                               m_probLo,
                               a_finestLevel,
                               a_ebGhost,
                               a_numGhost,
                               a_lsfGhost,
                               a_redistRad,
                               a_mgInterpOrder,
                               a_mgInterpRadius,
                               a_mgInterpWeight,
                               a_centroidInterpType,
                               a_ebInterpType,
                               m_baseif.at(phase::gas),
                               ebis_gas);

  m_realms[phase::solid]->define(a_grids,
                                 a_domains,
                                 a_refRat,
                                 a_dx,
                                 m_probLo,
                                 a_finestLevel,
                                 a_ebGhost,
                                 a_numGhost,
                                 a_lsfGhost,
                                 a_redistRad,
                                 a_mgInterpOrder,
                                 a_mgInterpRadius,
                                 a_mgInterpWeight,
                                 a_centroidInterpType,
                                 a_ebInterpType,
                                 m_baseif.at(phase::solid),
                                 ebis_sol);

  m_isDefined = true;
}

void
Realm::setGrids(const Vector<DisjointBoxLayout>& a_grids, const int a_finestLevel)
{
  CH_TIME("Realm::setGrids");
  if (m_verbosity > 5) {
    pout() << "Realm::setGrids" << endl;
  }

  for (auto& r : m_realms) {
    r.second->setGrids(a_grids, a_finestLevel);
  }
}

void
Realm::preRegrid()
{
  CH_TIME("Realm::preRegrid");
  if (m_verbosity > 5) {
    pout() << "Realm::preRegrid" << endl;
  }

  m_grids.resize(0);
  m_mflg.resize(0);
  m_validCells.resize(0);

  for (auto& mask : m_masks) {
    mask.second.resize(0);
  }
  //  m_masks.clear();

  for (auto& r : m_realms) {
    r.second->preRegrid();
  }
}

void
Realm::regridBase(const int a_lmin)
{
  CH_TIME("Realm::regridBase");
  if (m_verbosity > 5) {
    pout() << "Realm::regridBase" << endl;
  }

  for (auto& r : m_realms) {
    r.second->regridBase(a_lmin);
  }

  this->defineMFLevelGrid(a_lmin);
  this->defineValidCells();
  this->defineLevelTiles();
  this->defineMergeBlockGrid();
  this->defineParticleGhostMasks();
}

void
Realm::regridOperators(const int a_lmin)
{
  CH_TIME("Realm::regridOperators");
  if (m_verbosity > 5) {
    pout() << "Realm::regridOperators" << endl;
  }

  for (auto& r : m_realms) {
    r.second->regridOperators(a_lmin);
  }

  this->defineMasks(a_lmin);
  this->definePetscGrid();
}

void
Realm::defineMasks(const int a_lmin)
{
  CH_TIME("Realm::defineMasks");
  if (m_verbosity > 5) {
    pout() << "Realm::defineMasks" << endl;
  }

#ifdef CH_USE_PETSC
  // This mask is used when defining the PETSc interface and are therefore always required.
  this->registerMask(s_outer_cf_region, 1);
  this->registerMask(s_inner_cf_region, 1);
#endif

  // Regrid all masks
  this->defineOuterHaloMask(a_lmin);
  this->defineInnerHaloMask(a_lmin);

  this->defineOuterCFMask(a_lmin);
  this->defineInnerCFMask(a_lmin);

  this->defineCFIVS(a_lmin);
}

void
Realm::defineMFLevelGrid(const int /*a_lmin*/)
{
  CH_TIME("Realm::defineMFLevelGrid");
  if (m_verbosity > 5) {
    pout() << "Realm::defineMFLevelGrid" << endl;
  }

  m_mflg.resize(1 + m_finestLevel);
  m_mflgCoFi.resize(1 + m_finestLevel);
  m_mflgFiCo.resize(1 + m_finestLevel);

  PhaseRealm& gas = this->getRealm(phase::gas);
  PhaseRealm& sol = this->getRealm(phase::solid);

  const RefCountedPtr<EBIndexSpace>& ebis_gas = gas.getEBIndexSpace();
  const RefCountedPtr<EBIndexSpace>& ebis_sol = sol.getEBIndexSpace();

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    Vector<EBLevelGrid> eblgs;
    Vector<EBLevelGrid> eblgsCoFi;
    Vector<EBLevelGrid> eblgsFiCo;

    // Define the basic grids.
    //
    // The returned PhaseRealm::getEBLevelGridCoFi contains fine grids coarsened to the vector entry index. I.e.,
    // getEBLevelGridCoFi[0] contains the coarsening of grids on level 1.
    //
    // Similarly, PhaseRaelm::getEBLevelGridFiCo contains the refined grids. I.e., getEBLevelGrid[1] contains
    // the refinement of the grids on level 0. We dereference the pointers when we pass into MFLevelGrid,
    // so have to put some safeguards in place before we define.
    if (!ebis_gas.isNull()) {
      eblgs.push_back(*(gas.getEBLevelGrid()[lvl]));

      if (lvl < m_finestLevel) {
        eblgsCoFi.push_back(*(gas.getEBLevelGridCoFi()[lvl]));
      }
      if (lvl > 0) {
        eblgsFiCo.push_back(*(gas.getEBLevelGridFiCo()[lvl]));
      }
    }
    if (!ebis_sol.isNull()) {
      eblgs.push_back(*(sol.getEBLevelGrid()[lvl]));
      if (lvl < m_finestLevel) {
        eblgsCoFi.push_back(*(sol.getEBLevelGridCoFi()[lvl]));
      }
      if (lvl > 0) {
        eblgsFiCo.push_back(*(sol.getEBLevelGridFiCo()[lvl]));
      }
    }

    // Actual grids
    m_mflg[lvl] = RefCountedPtr<MFLevelGrid>(new MFLevelGrid(m_multifluidIndexSpace, eblgs));

    // Coarsened grids -- there nothing on the finest level since there's no finer grid.
    if (lvl < m_finestLevel) {
      m_mflgCoFi[lvl] = RefCountedPtr<MFLevelGrid>(new MFLevelGrid(m_multifluidIndexSpace, eblgsCoFi));
    }
    else {
      m_mflgCoFi[lvl] = RefCountedPtr<MFLevelGrid>(nullptr);
    }

    // Refined grids -- there's nothing on level 0 because there was nothing to refined from.
    if (lvl > 0) {
      m_mflgFiCo[lvl] = RefCountedPtr<MFLevelGrid>(new MFLevelGrid(m_multifluidIndexSpace, eblgsFiCo));
    }
    else {
      m_mflgFiCo[lvl] = RefCountedPtr<MFLevelGrid>(nullptr);
    }
  }
}

void
Realm::defineOuterHaloMask(const int /*a_lmin*/)
{
  CH_TIME("Realm::defineOuterHaloMask");
  if (m_verbosity > 5) {
    pout() << "Realm::defineOuterHaloMask" << endl;
  }

  // Loop through all masks and do something about the halo masks only.
  for (auto& m : m_masks) {

    // Get mask identifier and buffer.
    const std::string which_mask = m.first.first;
    const int         buffer     = m.first.second;

    if (which_mask == s_outer_particle_halo) {
      if (buffer <= 0) {
        MayDay::Abort("Realm::defineOuterHaloMask -- cannot have buffer <= 0!");
      }

      AMRMask& mask = m.second;

      mask.resize(1 + m_finestLevel);

      for (int lvl = 0; lvl < m_finestLevel; lvl++) {
        const DisjointBoxLayout& gridsCoar = m_grids[lvl];
        const DisjointBoxLayout& gridsFine = m_grids[lvl + 1];

        const ProblemDomain& domainCoar = m_domains[lvl];
        const ProblemDomain& domainFine = m_domains[lvl + 1];

        const int ncomp = 1;

        mask[lvl] = RefCountedPtr<LevelData<BaseFab<bool>>>(
          new LevelData<BaseFab<bool>>(gridsCoar, ncomp, IntVect::Zero));

        this->defineOuterHaloMask(*mask[lvl],
                                  domainCoar,
                                  domainFine,
                                  gridsCoar,
                                  gridsFine,
                                  buffer,
                                  m_refinementRatios[lvl]);
      }

      // Must explicitly set this because we want to the finest level mask to be nullptr, but we could be removing a grid level and in that case
      // the old mask will remain in the vector after resizing. This fixes that.
      mask[m_finestLevel] = RefCountedPtr<LevelData<BaseFab<bool>>>(nullptr);
    }
  }
}

void
Realm::defineOuterHaloMask(LevelData<BaseFab<bool>>& a_coarMask,
                           const ProblemDomain&      a_domainCoar,
                           const ProblemDomain& /*a_domainFine*/,
                           const DisjointBoxLayout& a_gridsCoar,
                           const DisjointBoxLayout& a_gridsFine,
                           const int                a_buffer,
                           const int                a_refRat)
{
  CH_TIME("Realm::defineOuterHaloMask");
  if (m_verbosity > 5) {
    pout() << "Realm::defineOuterHaloMask" << endl;
  }

  // TLDR: This routine defines a "mask" of valid coarse-grid cells (i.e., not covered by a finer level) around a fine grid. The mask is a_buffer wide. It is
  //       created by first fetching the cells around on the fine grid. This is a local operation where we get all the cells around each patch and then subtract
  //       cells that overlap with other boxes. This set of cells is then put on a LevelData<FArrayBox> mask so we can copy the result to the coarse grid and
  //       set the mask.

  constexpr int comp  = 0;
  constexpr int ncomp = 1;

  const DataIterator& ditCoar  = a_gridsCoar.dataIterator();
  const int           nboxCoar = ditCoar.size();

  // First, reset the mask.
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nboxCoar; mybox++) {
    const DataIndex& din = ditCoar[mybox];

    a_coarMask[din].setVal(false);
  }

  // Ok, we need a particle halo.
  if (a_buffer > 0) {

    // Coarsen the fine grid and make a mask on the coarsened fine grid
    DisjointBoxLayout dblCoFi;
    coarsen(dblCoFi, a_gridsFine, a_refRat);

    IntVectSet halo;

    const DataIterator& ditCoFi  = dblCoFi.dataIterator();
    const int           nboxCoFi = ditCoFi.size();

    // Go through the cofi grid and set the halo to true
#pragma omp parallel for schedule(runtime) reduction(+ : halo)
    for (int mybox = 0; mybox < nboxCoFi; mybox++) {
      const DataIndex& din = ditCoFi[mybox];

      const Box coFiBox = dblCoFi[din];

      // Make IntVect set consisting of only ghost cells (a_buffer) for each box.
      Box grownBox = grow(coFiBox, a_buffer);
      grownBox &= a_domainCoar;

      // Subtract non-ghosted box.
      IntVectSet myHalo(grownBox);
      myHalo -= coFiBox;

      // Subtract non-ghosted neighboring boxes.
      NeighborIterator nit(dblCoFi); // Neighbor iterator
      for (nit.begin(din); nit.ok(); ++nit) {
        const Box neighborBox = dblCoFi[nit()];
        myHalo -= neighborBox;
      }

      // Add to halo.
      halo |= myHalo;
    }

    // TLDR: In the above, we found the coarse-grid cells surrounding the fine level, viewed from the fine grids.
    //       Below, we create that view from the coarse grid. We use a BoxLayoutData<FArrayBox> on the coarsened fine grid,
    //       whose "ghost cells" can added to the _actual_ coarse grid. We then loop through those cells and set the mask.
    LevelData<FArrayBox> coFiMask(dblCoFi, ncomp, a_buffer * IntVect::Unit);
    LevelData<FArrayBox> coarMask(a_gridsCoar, ncomp, IntVect::Zero);

    // Reset masks
#pragma omp parallel
    {
#pragma omp for schedule(runtime)
      for (int mybox = 0; mybox < nboxCoFi; mybox++) {
        const DataIndex& din = ditCoFi[mybox];

        coFiMask[din].setVal(0.0);
      }
#pragma omp for schedule(runtime)
      for (int mybox = 0; mybox < nboxCoar; mybox++) {
        const DataIndex& din = ditCoar[mybox];

        coarMask[din].setVal(0.0);
      }

      // Run through the halo and set the halo cells to 1 on the coarsened fine grids. Since dblCoFi was a coarsening of the fine
      // grid, the mask value in the valid region (i.e., not including ghosts) is always zero. There used to be a bug here because
      // we only iterated through that region, but obviously the halo masks will always be zero in that case...
#pragma omp for schedule(runtime)
      for (int mybox = 0; mybox < nboxCoFi; mybox++) {
        const DataIndex& din = ditCoFi[mybox];

        const Box        region  = coFiMask[din].box();
        const IntVectSet curHalo = halo & region;
        for (IVSIterator ivsit(curHalo); ivsit.ok(); ++ivsit) {
          coFiMask[din](ivsit(), comp) = 1.0;
        }
      }
    }

    // Add the result to the coarse grid.
    Copier copier;
    copier.ghostDefine(dblCoFi, a_gridsCoar, a_domainCoar, a_buffer * IntVect::Unit);
    const Interval interv(0, 0);
    coFiMask.copyTo(interv, coarMask, interv, copier, LDaddOp<FArrayBox>());

    // Run through the grids and make the boolean mask
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nboxCoar; mybox++) {
      const DataIndex& din = ditCoar[mybox];

      const Box box = a_gridsCoar[din];

      BaseFab<bool>&   boolMask = a_coarMask[din];
      const FArrayBox& realMask = coarMask[din];

      // Not auto-vectorizable: one-time (regrid) setup loop with a data-dependent conditional write
      // to a bool mask.
      auto kernel = [&](const IntVect& iv) -> void {
        if (realMask(iv, comp) > 0.0) {
          boolMask(iv, comp) = true;
        }
      };

      BoxLoops::loop<D_DECL(1, 1, 1)>(box, kernel);
    }
  }
}

void
Realm::defineInnerHaloMask(const int /*a_lmin*/)
{
  CH_TIME("Realm::defineInnerHaloMask");
  if (m_verbosity > 5) {
    pout() << "Realm::defineInnerHaloMask" << endl;
  }

  // Loop through all masks and do something about the halo masks only.
  for (auto& m : m_masks) {

    // Get mask identifier and buffer.
    const std::string which_mask = m.first.first;
    const int         buffer     = m.first.second;

    if (which_mask == s_inner_particle_halo) {
      if (buffer <= 0) {
        MayDay::Abort("Realm::defineInnerHaloMask -- cannot have buffer <= 0!");
      }

      AMRMask& amrMask = m.second;

      amrMask.resize(1 + m_finestLevel);

      const int comp    = 0;
      const int numComp = 1;

      for (int lvl = 0; lvl <= m_finestLevel; lvl++) {

        const DisjointBoxLayout& dbl      = m_grids[lvl];
        const DataIterator&      dit      = dbl.dataIterator();
        const int                numBoxes = dit.size();

        amrMask[lvl] = RefCountedPtr<LevelData<BaseFab<bool>>>(
          new LevelData<BaseFab<bool>>(dbl, numComp, IntVect::Zero));

#pragma omp parallel for schedule(runtime)
        for (int mybox = 0; mybox < numBoxes; mybox++) {
          const DataIndex& din = dit[mybox];

          (*amrMask[lvl])[din].setVal(false);
        }

        if (lvl > 0) {
          // TLDR: Below, we use the valid cells to figure out the region on the inside of the refinement boundary. In everything below, the "fine" grid
          //       is the current grid level, and the coarse grid is the grid level below.

          const int refToCoar = m_refinementRatios[lvl - 1];

          const LevelData<BaseFab<bool>>& validCellsCoar = (*m_validCells[lvl - 1]);

          const DisjointBoxLayout& gridsCoar = m_grids[lvl - 1];
          const DisjointBoxLayout& gridsFine = m_grids[lvl];

          const ProblemDomain& domainCoar = m_domains[lvl - 1];

          const DataIterator& ditCoar = gridsCoar.dataIterator();
          const DataIterator& ditFine = gridsFine.dataIterator();

          const int numBoxesCoar = ditCoar.size();

          Vector<Box> boxesCoar = gridsCoar.boxArray();
          Vector<Box> boxesFine = gridsFine.boxArray();

          const Vector<int> procsCoar = gridsCoar.procIDs();
          const Vector<int> procsFine = gridsFine.procIDs();

          // Grow all the coarse boxes by the buffer and create a box layout on the coarse level using these boxes.
          for (int i = 0; i < boxesCoar.size(); i++) {
            boxesCoar[i] = grow(boxesCoar[i], buffer) & domainCoar;
          }

          // Coarsen all the fine boxes and create a boxlayout on the coarse level using the fine-grid decomposition
          for (int i = 0; i < boxesFine.size(); i++) {
            boxesFine[i] = coarsen(boxesFine[i], refToCoar);
          }

          const BoxLayout boxLayoutCoar(boxesCoar, procsCoar);
          const BoxLayout boxLayoutCoFi(boxesFine, procsFine);

          BoxLayoutData<FArrayBox> coarLevelMask(boxLayoutCoar, 1);
          BoxLayoutData<FArrayBox> coFiLevelMask(boxLayoutCoFi, 1);

          // Set the coarse mask to false everywhere on the coarse level, except on the region just inside the refinement boundary.
#pragma omp parallel for schedule(runtime)
          for (int mybox = 0; mybox < numBoxesCoar; mybox++) {
            const DataIndex& dinCoar = ditCoar[mybox];
            const Box&       cellBox = gridsCoar[dinCoar];

            FArrayBox&           coarMask   = coarLevelMask[dinCoar];
            const BaseFab<bool>& validCells = validCellsCoar[dinCoar];

            coarMask.setVal(0.0);

            // Not auto-vectorizable: one-time (regrid) setup loop with a data-dependent branch plus an
            // inner box-grow scatter (flagging all cells within the buffer).
            auto flagGrownRegion = [&](const IntVect& iv) -> void {
              if (validCells(iv, comp)) {
                const Box grownBox = grow(Box(iv, iv), buffer) & domainCoar;

#if CH_SPACEDIM == 3
                for (int k = grownBox.smallEnd(2); k <= grownBox.bigEnd(2); k++) {
#endif
                  for (int j = grownBox.smallEnd(1); j <= grownBox.bigEnd(1); j++) {
                    for (int i = grownBox.smallEnd(0); i <= grownBox.bigEnd(0); i++) {
                      coarMask(IntVect(D_DECL(i, j, k)), comp) = 1.0;
                    }
                  }
#if CH_SPACEDIM == 3
                }
#endif
              }
            };

            BoxLoops::loop<D_DECL(1, 1, 1)>(cellBox, flagGrownRegion);
          }

          // Reset the mask on the coarsened grid
#pragma omp parallel for schedule(runtime)
          for (int mybox = 0; mybox < numBoxes; mybox++) {
            const DataIndex& din = dit[mybox];

            FArrayBox& coFiMask = coFiLevelMask[din];

            coFiMask.setVal(0.0);
          }

          // Increment the data from the coarse grid to the coarsened fine grid. This must
          // also increment with the ghost well values.
          const Interval srcInterv(comp, comp);
          const Interval dstInterv(comp, comp);

          coarLevelMask.addTo(srcInterv, coFiLevelMask, dstInterv, domainCoar);

          // Iterate through the current level and, for each cell, check if the coarsened cell value was
          // flagged as a mask cell.
#pragma omp parallel for schedule(runtime)
          for (int mybox = 0; mybox < numBoxes; mybox++) {
            const DataIndex& din     = dit[mybox];
            const Box&       cellBox = gridsFine[din];

            BaseFab<bool>&   curMask  = (*amrMask[lvl])[din];
            const FArrayBox& coFiMask = coFiLevelMask[din];

            // Not auto-vectorizable: one-time (regrid) setup loop — coarsen() gather plus a
            // conditional bool write.
            auto kernel = [&](const IntVect& iv) -> void {
              const IntVect coarIV = coarsen(iv, refToCoar);

              if (coFiMask(coarIV, comp) > 0.0) {
                curMask(iv, comp) = true;
              }
            };

            BoxLoops::loop<D_DECL(1, 1, 1)>(cellBox, kernel);
          }
        }
      }
    }
  }
}

void
Realm::defineOuterCFMask(const int /*a_lmin*/)
{
  CH_TIME("Realm::defineOuterCFMask");
  if (m_verbosity > 5) {
    pout() << "Realm::defineOuterCFMask" << endl;
  }

  // Loop through all masks and do something about the halo masks only.
  for (auto& m : m_masks) {

    // Get mask identifier and buffer.
    const std::string which_mask = m.first.first;
    const int         buffer     = m.first.second;

    if (which_mask == s_outer_cf_region) {
      if (buffer <= 0) {
        MayDay::Abort("Realm::defineOuterCFMask -- cannot have buffer <= 0!");
      }

      AMRMask& mask = m.second;

      mask.resize(1 + m_finestLevel);

      for (int lvl = 0; lvl < m_finestLevel; lvl++) {
        const DisjointBoxLayout& gridsCoar = m_grids[lvl];
        const DisjointBoxLayout& gridsFine = m_grids[lvl + 1];

        const ProblemDomain& domainCoar = m_domains[lvl];
        const ProblemDomain& domainFine = m_domains[lvl + 1];

        const int ncomp = 1;

        mask[lvl] = RefCountedPtr<LevelData<BaseFab<bool>>>(
          new LevelData<BaseFab<bool>>(gridsCoar, ncomp, IntVect::Zero));

        this->defineOuterCFMask(*mask[lvl],
                                domainCoar,
                                domainFine,
                                gridsCoar,
                                gridsFine,
                                buffer,
                                m_refinementRatios[lvl]);
      }

      // Must explicitly set this because we want to the finest level mask to be nullptr, but we could be removing a grid level and in that case
      // the old mask will remain in the vector after resizing. This fixes that.
      mask[m_finestLevel] = RefCountedPtr<LevelData<BaseFab<bool>>>(nullptr);
    }
  }
}

void
Realm::defineOuterCFMask(LevelData<BaseFab<bool>>& a_coarMask,
                         const ProblemDomain&      a_domainCoar,
                         const ProblemDomain& /*a_domainFine*/,
                         const DisjointBoxLayout& a_gridsCoar,
                         const DisjointBoxLayout& a_gridsFine,
                         const int                a_buffer,
                         const int                a_refRat)
{
  CH_TIME("Realm::defineOuterCFMask");
  if (m_verbosity > 5) {
    pout() << "Realm::defineOuterCFMask" << endl;
  }

  // TLDR: This routine defines a "mask" of valid coarse-grid cells (i.e., not covered by a finer level) around a fine grid. The mask is a_buffer wide. It is
  //       created by first fetching the cells around on the fine grid. This is a local operation where we get all the cells around each patch and then subtract
  //       cells that overlap with other boxes. This set of cells is then put on a LevelData<FArrayBox> mask so we can copy the result to the coarse grid and
  //       set the mask.

  constexpr int comp  = 0;
  constexpr int ncomp = 1;

  const DataIterator& ditCoar  = a_gridsCoar.dataIterator();
  const int           nboxCoar = ditCoar.size();

  // First, reset the mask.
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nboxCoar; mybox++) {
    const DataIndex& din = ditCoar[mybox];

    a_coarMask[din].setVal(false);
  }

  // Ok, we need a particle halo.
  if (a_buffer > 0) {

    // Coarsen the fine grid and make a mask on the coarsened fine grid
    DisjointBoxLayout dblCoFi;
    coarsen(dblCoFi, a_gridsFine, a_refRat);

    IntVectSet halo;

    const DataIterator& ditCoFi  = dblCoFi.dataIterator();
    const int           nboxCoFi = ditCoFi.size();

    // Go through the cofi grid and set the halo to true
#pragma omp parallel for schedule(runtime) reduction(+ : halo)
    for (int mybox = 0; mybox < nboxCoFi; mybox++) {
      const DataIndex& din = ditCoFi[mybox];

      const Box coFiBox = dblCoFi[din];

      // Make IntVect set consisting of only ghost cells (a_buffer) for each box.
      Box grownBox = grow(coFiBox, a_buffer);
      grownBox &= a_domainCoar;

      // Subtract non-ghosted box.
      IntVectSet myHalo(grownBox);
      myHalo -= coFiBox;

      // Subtract non-ghosted neighboring boxes.
      NeighborIterator nit(dblCoFi); // Neighbor iterator
      for (nit.begin(din); nit.ok(); ++nit) {
        const Box neighborBox = dblCoFi[nit()];
        myHalo -= neighborBox;
      }

      // Add to halo.
      halo |= myHalo;
    }

    // TLDR: In the above, we found the coarse-grid cells surrounding the fine level, viewed from the fine grids.
    //       Below, we create that view from the coarse grid. We use a BoxLayoutData<FArrayBox> on the coarsened fine grid,
    //       whose "ghost cells" can added to the _actual_ coarse grid. We then loop through those cells and set the mask.
    LevelData<FArrayBox> coFiMask(dblCoFi, ncomp, a_buffer * IntVect::Unit);
    LevelData<FArrayBox> coarMask(a_gridsCoar, ncomp, IntVect::Zero);

    // Reset masks
#pragma omp parallel
    {
#pragma omp for schedule(runtime)
      for (int mybox = 0; mybox < nboxCoFi; mybox++) {
        const DataIndex& din = ditCoFi[mybox];

        coFiMask[din].setVal(0.0);
      }
#pragma omp for schedule(runtime)
      for (int mybox = 0; mybox < nboxCoar; mybox++) {
        const DataIndex& din = ditCoar[mybox];

        coarMask[din].setVal(0.0);
      }

      // Run through the halo and set the halo cells to 1 on the coarsened fine grids. Since dblCoFi was a coarsening of the fine
      // grid, the mask value in the valid region (i.e., not including ghosts) is always zero. There used to be a bug here because
      // we only iterated through that region, but obviously the halo masks will always be zero in that case...
#pragma omp for schedule(runtime)
      for (int mybox = 0; mybox < nboxCoFi; mybox++) {
        const DataIndex& din = ditCoFi[mybox];

        const Box        region  = coFiMask[din].box();
        const IntVectSet curHalo = halo & region;
        for (IVSIterator ivsit(curHalo); ivsit.ok(); ++ivsit) {
          coFiMask[din](ivsit(), comp) = 1.0;
        }
      }
    }

    // Add the result to the coarse grid.
    Copier copier;
    copier.ghostDefine(dblCoFi, a_gridsCoar, a_domainCoar, a_buffer * IntVect::Unit);
    const Interval interv(0, 0);
    coFiMask.copyTo(interv, coarMask, interv, copier, LDaddOp<FArrayBox>());

    // Run through the grids and make the boolean mask
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nboxCoar; mybox++) {
      const DataIndex& din = ditCoar[mybox];

      const Box box = a_gridsCoar[din];

      BaseFab<bool>&   boolMask = a_coarMask[din];
      const FArrayBox& realMask = coarMask[din];

      // Not auto-vectorizable: one-time (regrid) setup loop with a data-dependent conditional write
      // to a bool mask.
      auto kernel = [&](const IntVect& iv) -> void {
        if (realMask(iv, comp) > 0.0) {
          boolMask(iv, comp) = true;
        }
      };

      BoxLoops::loop<D_DECL(1, 1, 1)>(box, kernel);
    }
  }
}

void
Realm::defineInnerCFMask(const int /*a_lmin*/)
{
  CH_TIME("Realm::defineInnerCFMask");
  if (m_verbosity > 5) {
    pout() << "Realm::defineInnerCFMask" << endl;
  }

  const int numComp = 1;

  // Loop through all masks and do something about the halo masks only.
  for (auto& m : m_masks) {
    const std::string which_mask = m.first.first;
    const int         buffer     = m.first.second;

    if (which_mask == s_inner_cf_region) {
      if (buffer <= 0) {
        MayDay::Abort("Realm::defineInnerCFMask -- cannot have buffer <= 0!");
      }

      AMRMask& amrMask = m.second;

      amrMask.resize(1 + m_finestLevel);

      for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
        const DisjointBoxLayout& dbl      = m_grids[lvl];
        const ProblemDomain&     domain   = m_domains[lvl];
        const DataIterator&      dit      = dbl.dataIterator();
        const int                numBoxes = dit.size();

        amrMask[lvl] = RefCountedPtr<LevelData<BaseFab<bool>>>(
          new LevelData<BaseFab<bool>>(dbl, numComp, IntVect::Zero));

#pragma omp parallel for schedule(runtime)
        for (int mybox = 0; mybox < numBoxes; mybox++) {
          const DataIndex& din     = dit[mybox];
          const Box        cellBox = dbl[din];

          BaseFab<bool>& mask = (*amrMask[lvl])[din];

          mask.setVal(false);

          if (lvl > 0) {
            Vector<Box> neighborBoxes(1, cellBox);

            NeighborIterator nit(dbl);

            for (nit.begin(din); nit.ok(); ++nit) {
              neighborBoxes.push_back(dbl[nit()]);
            }

            for (int dir = 0; dir < SpaceDim; dir++) {
              const Box sideBoxLo = adjCellLo(cellBox, dir, -buffer);
              const Box sideBoxHi = adjCellHi(cellBox, dir, -buffer);

              // Not auto-vectorizable: one-time (regrid) setup loop — per-cell box-grow and
              // neighbor-overlap counting (box arithmetic plus an inner loop over neighbor boxes).
              auto kernel = [&](const IntVect& iv) -> void {
                const Box box = grow(Box(iv, iv), buffer) & domain;

                int c = static_cast<int>(box.numPts());

                for (int i = 0; i < neighborBoxes.size(); i++) {
                  const Box overlapBox = box & neighborBoxes[i];

                  c -= static_cast<int>(overlapBox.numPts());
                }

                mask(iv) = (c == 0) ? false : true;
              };

              BoxLoops::loop<D_DECL(1, 1, 1)>(sideBoxLo, kernel);
              BoxLoops::loop<D_DECL(1, 1, 1)>(sideBoxHi, kernel);
            }
          }
        }
      }
    }
  }
}

void
Realm::defineCFIVS(const int /*a_lmin*/)
{
  CH_TIME("Realm::defineCFIVS");
  if (m_verbosity > 5) {
    pout() << "Realm::defineCFIVS" << endl;
  }

  // Loop through all masks and do something about the halo masks only.
  for (auto& m : m_masks) {

    // Get mask identifier and buffer.
    const std::string which_mask = m.first.first;
    const int         buffer     = m.first.second;

    if (which_mask == s_cfivs) {
      if (buffer < 0) {
        MayDay::Abort("Realm::defineCFIVS -- cannot have buffer <= 0!");
      }

      AMRMask& amrMask = m.second;

      amrMask.resize(1 + m_finestLevel);

      const int comp    = 0;
      const int numComp = 1;

      for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
        const ProblemDomain&     domain   = m_domains[lvl];
        const DisjointBoxLayout& dbl      = m_grids[lvl];
        const DataIterator&      dit      = dbl.dataIterator();
        const int                numBoxes = dit.size();

        amrMask[lvl] = RefCountedPtr<LevelData<BaseFab<bool>>>(
          new LevelData<BaseFab<bool>>(dbl, numComp, m_numGhost * IntVect::Unit));

#pragma omp parallel for schedule(runtime)
        for (int mybox = 0; mybox < numBoxes; mybox++) {
          const DataIndex& din      = dit[mybox];
          const Box        cellBox  = dbl[din];
          const Box        ghostBox = grow(cellBox, m_numGhost) & domain;

          BaseFab<bool>& mask = (*amrMask[lvl])[din];

          if (lvl == 0) {
            mask.setVal(false);
          }
          else {

            mask.setVal(true, ghostBox, comp);
            mask.setVal(false, cellBox, comp);

            NeighborIterator nit(dbl);
            for (nit.begin(din); nit.ok(); ++nit) {
              const Box neighborBox = dbl[nit()];

              mask.setVal(false, neighborBox, comp);
            }
          }
        }
      }
    }
  }
}

void
Realm::defineValidCells()
{
  CH_TIME("Realm::defineValidCells");
  if (m_verbosity > 5) {
    pout() << "Realm::defineValidCells" << endl;
  }

  constexpr int curComp  = 0;
  constexpr int numComp  = 1;
  constexpr int numGhost = 0;

  m_validCells.resize(1 + m_finestLevel);

  for (int lvl = m_finestLevel; lvl >= 0; lvl--) {
    const DisjointBoxLayout& dblCoar  = m_grids[lvl];
    const DataIterator&      ditCoar  = dblCoar.dataIterator();
    const int                nboxCoar = ditCoar.size();

    m_validCells[lvl] = RefCountedPtr<LevelData<BaseFab<bool>>>(
      new LevelData<BaseFab<bool>>(dblCoar, numComp, numGhost * IntVect::Unit));

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nboxCoar; mybox++) {
      const DataIndex& din = ditCoar[mybox];

      BaseFab<bool>& validCells = (*m_validCells[lvl])[din];

      validCells.setVal(true);
    }

    if (lvl < m_finestLevel) {
      // If there is a finer level then we need to coarsen that level onto the current level and set all those cells to
      // false. I'm calling the current level the 'coarse level' and the finer level the 'fine level'.

      // Coarsened fine grids.
      DisjointBoxLayout dblCoFi;
      coarsen(dblCoFi, m_grids[lvl + 1], m_refinementRatios[lvl]);

      // Create some data = 1 on the fine grid and = 0 on the coarse grid
      LevelData<FArrayBox> coFiData(dblCoFi, numComp, numGhost * IntVect::Unit);
      LevelData<FArrayBox> coarData(dblCoar, numComp, numGhost * IntVect::Unit);

      const DataIterator& ditCoFi  = dblCoFi.dataIterator();
      const int           nboxCoFi = ditCoFi.size();

#pragma omp parallel
      {
#pragma omp for schedule(runtime)
        for (int mybox = 0; mybox < nboxCoFi; mybox++) {
          const DataIndex& din = ditCoFi[mybox];

          coFiData[din].setVal(1.0);
        }

#pragma omp for schedule(runtime)
        for (int mybox = 0; mybox < nboxCoar; mybox++) {
          const DataIndex& din = ditCoar[mybox];

          coarData[din].setVal(0.0);
        }
      }

      // Copy from fine to coarse.
      Copier copier(dblCoFi, dblCoar);
      coFiData.copyTo(Interval(0, 0), coarData, Interval(0, 0), copier);

      // Go through the coarse grid -- wherever we find a value of 1 there is also a fine grid.
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nboxCoar; mybox++) {
        const DataIndex& din = ditCoar[mybox];

        BaseFab<bool>&   boolMask = (*m_validCells[lvl])[din];
        const FArrayBox& fabMask  = coarData[din];

        // Not auto-vectorizable: one-time (regrid) setup loop with a data-dependent conditional write
        // to a bool mask.
        auto kernel = [&](const IntVect& iv) -> void {
          if (fabMask(iv, curComp) > 0.0) {
            boolMask(iv, curComp) = false;
          }
        };

        BoxLoops::loop<D_DECL(1, 1, 1)>(fabMask.box(), kernel);
      }
    }
  }
}

void
Realm::registerParticleGhostMask(const int a_width) noexcept
{
  CH_TIME("Realm::registerParticleGhostMask");
  if (m_verbosity > 5) {
    pout() << "Realm::registerParticleGhostMask" << endl;
  }

  if (a_width <= 0) {
    MayDay::Abort("Realm::registerParticleGhostMask -- particle ghost-mask width must be > 0");
  }

  m_particleGhostMaskWidths.insert(a_width);
}

void
Realm::defineParticleGhostMasks() noexcept
{
  CH_TIME("Realm::defineParticleGhostMasks");
  if (m_verbosity > 5) {
    pout() << "Realm::defineParticleGhostMasks" << endl;
  }

  m_particleGhostMask.clear();
  m_particleGhostMaskFineToCoar.clear();
  m_particleGhostMaskCoarToFine.clear();

  // Particle ghost filling is intentionally NOT supported on periodic domains (a periodic ghost would
  // require wrapping particle copies across the domain). Abort rather than silently mis-fill if a width
  // is registered on a periodic domain.
  if (!m_particleGhostMaskWidths.empty()) {
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      if (m_domains[lvl].isPeriodic()) {
        MayDay::Abort("Realm::defineParticleGhostMasks -- particle ghost masks do not support periodic domains");
      }
    }
  }

  // No particle ghost masks are built unless downstream code has registered at least one width.
  for (const int ghost : m_particleGhostMaskWidths) {
    if (ghost >= m_minBlockSize) {
      const std::string msg = "Realm::defineParticleGhostMasks -- ghost width = " + std::to_string(ghost) +
                              " must be < min_block_size = " + std::to_string(m_minBlockSize);
      MayDay::Abort(msg.c_str());
    }

    AMRParticleGhostMask& same       = m_particleGhostMask[ghost];
    AMRParticleGhostMask& fineToCoar = m_particleGhostMaskFineToCoar[ghost];
    AMRParticleGhostMask& coarToFine = m_particleGhostMaskCoarToFine[ghost];

    same.resize(1 + m_finestLevel);
    fineToCoar.resize(1 + m_finestLevel);
    coarToFine.resize(1 + m_finestLevel);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const DisjointBoxLayout& dbl    = m_grids[lvl];
      const ProblemDomain&     domain = m_domains[lvl];

      same[lvl]       = RefCountedPtr<LayoutData<ParticleGhostMask>>(new LayoutData<ParticleGhostMask>(dbl));
      fineToCoar[lvl] = RefCountedPtr<LayoutData<ParticleGhostMask>>(new LayoutData<ParticleGhostMask>(dbl));
      coarToFine[lvl] = RefCountedPtr<LayoutData<ParticleGhostMask>>(new LayoutData<ParticleGhostMask>(dbl));

      // SAME level: targets are the abutting boxes on this level (each grown by 'ghost' at this resolution).
      this->defineParticleGhostMaskSameLevel(*same[lvl], dbl, domain, ghost);

      // COARSER (lvl -> lvl-1): fine cells that scatter DOWN to the coarse level. Restricted to the fine
      // side of the coarse-fine halo; ghost width in destination (coarse) cells. CONTRACT: fineToCoar[0]
      // is intentionally left allocated-but-undefined (level 0 has no coarser level); callers must not
      // query it (gatherGhostsFromMasks guards with a_srcLvl >= 1).
      if (lvl > 0) {
        this->defineParticleGhostMaskFineToCoar(*fineToCoar[lvl],
                                                dbl,
                                                m_grids[lvl - 1],
                                                m_domains[lvl - 1],
                                                m_refinementRatios[lvl - 1],
                                                ghost);
      }

      // FINER (lvl -> lvl+1): coarse cells that scatter UP to the finer level. Restricted to the coarse
      // side of the coarse-fine halo; ghost width in destination (fine) cells. CONTRACT: coarToFine
      // [m_finestLevel] is intentionally left allocated-but-undefined (the finest level has no finer
      // level); callers must not query it (gatherGhostsFromMasks guards with a_srcLvl < m_finestLevel).
      if (lvl < m_finestLevel) {
        this->defineParticleGhostMaskCoarToFine(*coarToFine[lvl],
                                                dbl,
                                                m_grids[lvl + 1],
                                                *m_validCells[lvl],
                                                m_domains[lvl + 1],
                                                m_refinementRatios[lvl],
                                                ghost);
      }
    }
  }
}

void
Realm::defineParticleGhostMaskSameLevel(LayoutData<ParticleGhostMask>& a_mask,
                                        const DisjointBoxLayout&       a_dbl,
                                        const ProblemDomain&           a_domain,
                                        const int                      a_ghost) noexcept
{
  CH_TIME("Realm::defineParticleGhostMaskSameLevel");
  if (m_verbosity > 5) {
    pout() << "Realm::defineParticleGhostMaskSameLevel" << endl;
  }

  // A same-level neighbour N contributes a target to this box's cells that lie within a_ghost of N, i.e.
  // grow(N, a_ghost) & box. NeighborIterator yields the abutting boxes (never self), so the whole build
  // is box-local: each box fills its own CSR with no cross-box writes. Two passes per box (count, pack).
  // Thread-safe over boxes: distinct boxes write disjoint LayoutData storage, all neighbour lookups are
  // read-only, and NeighborIterator (mutable cursor) is declared per iteration so it is thread-private.
  const DataIterator& dit  = a_dbl.dataIterator();
  const int           nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex&   din = dit[mybox];
    const Box          box = a_dbl[din];
    ParticleGhostMask& t   = a_mask[din];

    NeighborIterator nit(a_dbl);

    t.define(box);

    // Count targets per cell.
    for (nit.begin(din); nit.ok(); ++nit) {
      const Box reach = grow(a_dbl[nit()], a_ghost) & box;
      for (BoxIterator bit(reach); bit.ok(); ++bit) {
        t.incrementCount(bit());
      }
    }

    t.allocate();

    // Pack the targets.
    for (nit.begin(din); nit.ok(); ++nit) {
      const LevelTiles::BoxIDs target(a_dbl.index(nit()), a_dbl.procID(nit()));
      const Box                reach = grow(a_dbl[nit()], a_ghost) & box;

      for (BoxIterator bit(reach); bit.ok(); ++bit) {
        t.addTarget(bit(), target);
      }
    }

    t.finalize();
  }
}

void
Realm::defineParticleGhostMaskFineToCoar(LayoutData<ParticleGhostMask>& a_mask,
                                         const DisjointBoxLayout&       a_thisLayout,
                                         const DisjointBoxLayout&       a_coarLayout,
                                         const ProblemDomain&           a_coarDomain,
                                         const int                      a_refRat,
                                         const int                      a_ghost) noexcept
{
  CH_TIME("Realm::defineParticleGhostMaskFineToCoar");
  if (m_verbosity > 5) {
    pout() << "Realm::defineParticleGhostMaskFineToCoar" << endl;
  }

  // Work at the COARSE (target) resolution by coarsening this (fine) level; the ghost width is then
  // measured in coarse cells. coarsen() preserves this level's DataIndex, so a motion item's toIndex also
  // indexes a_mask/a_thisLayout.
  DisjointBoxLayout coFiLayout;
  coarsen(coFiLayout, a_thisLayout, a_refRat);

  // ghostDefine: each coarse box's grown region (a_ghost coarse cells) -> the coarsened-fine footprint.
  // Gives, per footprint cell, ALL coarse boxes whose ghosted box contains it (multi-target).
  Copier copier;
  copier.ghostDefine(a_coarLayout, coFiLayout, a_coarDomain, a_ghost * IntVect::Unit);

  const DataIterator& dit  = a_thisLayout.dataIterator();
  const int           nbox = dit.size();

  // Restrict sources to the INNER coarse-fine halo: footprint cells within a_ghost coarse cells of the
  // coarse-fine interface, built geometrically from the coarsened-fine footprint. Stored per box as a
  // DenseIntVectSet (a bitmask over the footprint box), so the build needs no BaseFab conversion.
  //
  // This loop does both per-box setups at once (init the box's CSR a_mask[din] and its halo innerHalo[din]).
  // Thread-safe over boxes: distinct boxes write disjoint storage, all layout lookups are read-only, and
  // the DenseIntVectSet scratch and the NeighborIterator are declared per iteration (thread-private).
  // coFiLayout shares this level's DataIndex (coarsen preserves it), so dit/nbox iterate it too.
  LayoutData<DenseIntVectSet> innerHalo(coFiLayout);

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din      = dit[mybox];
    const Box        coFiBox  = coFiLayout[din];
    const Box        grownBox = grow(coFiBox, a_ghost) & a_coarDomain;

    a_mask[din].define(a_thisLayout[din]);

    // Ring = the VALID coarse cells surrounding this footprint box: reach OUT to grownBox, then remove the
    // footprint box itself and its abutting footprint neighbours (internal seams).
    DenseIntVectSet ring(grownBox, true);
    ring -= coFiBox;

    NeighborIterator nit(coFiLayout);
    for (nit.begin(din); nit.ok(); ++nit) {
      ring -= coFiLayout[nit()];
    }

    // Inner halo = footprint cells within a_ghost of that valid ring: grow the ring back IN, then keep the
    // footprint cells (intersect a full DenseIntVectSet over coFiBox, which also fixes the halo's domain to
    // the footprint box). The ring MUST be the subtracted one -- otherwise every footprint cell is kept.
    ring.grow(a_ghost);

    DenseIntVectSet& halo = innerHalo[din];
    halo                  = DenseIntVectSet(coFiBox, true);
    halo &= ring;
  }

  // Each retained (inner-halo) footprint cell refines to whole fine source cells; every such cell keeps
  // ALL coarse target boxes reaching it. Two motion-plan walks (count, then packed fill). These stay
  // serial: distinct motion items can target the same destination box, so incrementCount/addTarget would
  // race, and they iterate the Copier plan rather than the grid boxes.
  const auto forEachContribution = [&](auto&& a_emit) {
    const CopyIterator::local_from_to plans[2] = {CopyIterator::LOCAL, CopyIterator::TO};

    for (const auto plan : plans) {
      for (CopyIterator cit(copier, plan); cit.ok(); ++cit) {
        const MotionItem&        item = cit();
        const LevelTiles::BoxIDs target(a_coarLayout.index(item.fromIndex), a_coarLayout.procID(item.fromIndex));
        const Box                thisBox = a_thisLayout[item.toIndex]; // fine source box
        const DenseIntVectSet&   halo    = innerHalo[item.toIndex];

        for (BoxIterator bit(item.toRegion); bit.ok(); ++bit) {
          const IntVect cc = bit();
          if (!halo[cc]) {
            continue;
          }

          const Box fineCells = refine(Box(cc, cc), a_refRat) & thisBox;

          for (BoxIterator fit(fineCells); fit.ok(); ++fit) {
            a_emit(item.toIndex, fit(), target);
          }
        }
      }
    }
  };

  forEachContribution([&](const DataIndex& a_di, const IntVect& a_iv, const LevelTiles::BoxIDs&) {
    a_mask[a_di].incrementCount(a_iv);
  });

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    a_mask[dit[mybox]].allocate();
  }

  forEachContribution([&](const DataIndex& a_di, const IntVect& a_iv, const LevelTiles::BoxIDs& a_tg) {
    a_mask[a_di].addTarget(a_iv, a_tg);
  });

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    a_mask[dit[mybox]].finalize();
  }
}

void
Realm::defineParticleGhostMaskCoarToFine(LayoutData<ParticleGhostMask>&  a_mask,
                                         const DisjointBoxLayout&        a_thisLayout,
                                         const DisjointBoxLayout&        a_fineLayout,
                                         const LevelData<BaseFab<bool>>& a_validCells,
                                         const ProblemDomain&            a_fineDomain,
                                         const int                       a_refRat,
                                         const int                       a_ghost) noexcept
{
  CH_TIME("Realm::defineParticleGhostMaskCoarToFine");
  if (m_verbosity > 5) {
    pout() << "Realm::defineParticleGhostMaskCoarToFine" << endl;
  }

  // Work at the FINE (target) resolution by refining this (coarse) level; the ghost width is then measured
  // in fine cells. refine() preserves this level's DataIndex, so a motion item's toIndex also indexes
  // a_mask/a_thisLayout.
  DisjointBoxLayout refinedLayout;
  refine(refinedLayout, a_thisLayout, a_refRat);

  // ghostDefine: each fine box's grown region (a_ghost fine cells) -> this level's cells at fine res. Gives,
  // per fine cell, ALL fine boxes whose ghosted box contains it (multi-target).
  Copier copier;
  copier.ghostDefine(a_fineLayout, refinedLayout, a_fineDomain, a_ghost * IntVect::Unit);

  const DataIterator& dit  = a_thisLayout.dataIterator();
  const int           nbox = dit.size();

  // Per-box CSR setup/teardown loops are embarrassingly parallel (each writes only its own box). The two
  // motion-plan walks below stay serial: distinct motion items can target the same destination box, so
  // their incrementCount/addTarget would race, and they iterate the Copier plan rather than the grid boxes.
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    a_mask[dit[mybox]].define(a_thisLayout[dit[mybox]]);
  }

  // Coarsen the fine reach to whole coarse source cells and keep only VALID ones (covered cells have no
  // particles) -- that is the coarse-side outer halo; every retained cell keeps ALL fine target boxes
  // reaching it. Each contribution also carries an ACCEPTANCE BOX: the fine box's ghost region (grown by
  // a_ghost fine cells). The scatter uses it to prune whole-coarse-cell over-communication down to the
  // exact fine shell (a coarse particle is scattered only if its FINE cell lies in the acceptance box).
  // Domains are non-periodic here (see defineParticleGhostMasks), so the fine and this-level index spaces
  // coincide -- no periodic shift is needed. Two motion-plan walks.
  const auto forEachContribution = [&](auto&& a_emit) {
    const CopyIterator::local_from_to plans[2] = {CopyIterator::LOCAL, CopyIterator::TO};

    for (const auto plan : plans) {
      for (CopyIterator cit(copier, plan); cit.ok(); ++cit) {
        const MotionItem&        item = cit();
        const LevelTiles::BoxIDs target(a_fineLayout.index(item.fromIndex), a_fineLayout.procID(item.fromIndex));
        const Box                thisBox     = a_thisLayout[item.toIndex]; // coarse source box
        const BaseFab<bool>&     valid       = a_validCells[item.toIndex];
        const Box                coarseCells = coarsen(item.toRegion, a_refRat) & thisBox;

        // Acceptance box in fine (destination) cells: the fine box grown by the ghost width.
        const Box acceptBox = grow(a_fineLayout[item.fromIndex], a_ghost);

        for (BoxIterator bit(coarseCells); bit.ok(); ++bit) {
          const IntVect cc = bit();
          if (!valid(cc, 0)) {
            continue;
          }
          a_emit(item.toIndex, cc, target, acceptBox);
        }
      }
    }
  };

  forEachContribution([&](const DataIndex& a_di, const IntVect& a_iv, const LevelTiles::BoxIDs&, const Box&) {
    a_mask[a_di].incrementCount(a_iv);
  });

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    a_mask[dit[mybox]].allocate();
    a_mask[dit[mybox]].allocateTargetBoxes();
  }

  forEachContribution(
    [&](const DataIndex& a_di, const IntVect& a_iv, const LevelTiles::BoxIDs& a_tg, const Box& a_acceptBox) {
      a_mask[a_di].addTarget(a_iv, a_tg, a_acceptBox);
    });

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    a_mask[dit[mybox]].finalize();
  }
}

void
Realm::defineLevelTiles() noexcept
{
  CH_TIME("Realm::defineLevelTiles");
  if (m_verbosity > 5) {
    pout() << "Realm::defineLevelTiles" << endl;
  }

  m_levelTiles.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    m_levelTiles[lvl] = RefCountedPtr<LevelTiles>(new LevelTiles(m_grids[lvl], m_minBlockSize));
  }
}

void
Realm::defineMergeBlockGrid() noexcept
{
  CH_TIME("Realm::defineMergeBlockGrid");
  if (m_verbosity > 5) {
    pout() << "Realm::defineMergeBlockGrid" << endl;
  }

  m_mergeBlockGrids.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    m_mergeBlockGrids[lvl] = RefCountedPtr<MergeBlockGrid>(new MergeBlockGrid(m_grids[lvl], m_minBlockSize));
  }
}

void
Realm::definePetscGrid() noexcept
{
#ifdef CH_USE_PETSC
  CH_TIME("Realm::definePetscGrid");
  if (m_verbosity > 5) {
    pout() << "Realm::definePetscGrid" << endl;
  }

  if (!(m_petscGrid.isNull())) {
    m_petscGrid->clear();
  }
  else {
    m_petscGrid = RefCountedPtr<PetscGrid>(new PetscGrid());
  }

  m_petscGrid->define(m_mflg,
                      m_mflgCoFi,
                      m_mflgFiCo,
                      m_validCells,
                      this->getMask(s_outer_cf_region, 1),
                      this->getMask(s_inner_cf_region, 1),
                      m_refinementRatios,
                      m_dx,
                      m_finestLevel,
                      m_numGhost);
#endif
}

void
Realm::registerOperator(const std::string& a_operator, const phase::which_phase a_phase)
{
  CH_TIME("Realm::registerOperator(operator, phase)");
  if (m_verbosity > 5) {
    pout() << "Realm::registerOperator(operator, phase)" << endl;
  }

  m_realms[a_phase]->registerOperator(a_operator);
}

bool
Realm::queryOperator(const std::string& a_operator, const phase::which_phase a_phase) const
{
  CH_TIME("Realm::queryOperator");
  if (m_verbosity > 5) {
    pout() << "Realm::queryOperator" << endl;
  }

  return m_realms[a_phase]->queryOperator(a_operator);
}

void
Realm::registerMask(const std::string& a_mask, const int a_buffer)
{
  CH_TIME("Realm::registerMask(mask, buffer)");
  if (m_verbosity > 5) {
    pout() << "Realm::registerMask(mask, buffer)" << endl;
  }

  m_masks.emplace(std::pair<string, int>(a_mask, a_buffer), AMRMask());
}

bool
Realm::queryMask(const std::string& a_mask, const int a_buffer) const
{
  CH_TIME("Realm::queryMask(mask, buffer)");
  if (m_verbosity > 5) {
    pout() << "Realm::queryMask(mask, buffer)" << endl;
  }

  bool ret;

  if (m_masks.count(std::pair<string, int>(a_mask, a_buffer)) > 0) {
    ret = true;
  }
  else {
    ret = false;
  }

  return ret;
}

PhaseRealm&
Realm::getRealm(const phase::which_phase a_phase)
{
  return *m_realms[a_phase];
}

const Vector<int>&
Realm::getRefinementRatios() const
{
  return m_refinementRatios;
}

const Vector<Real>&
Realm::getDx() const
{
  return m_dx;
}

const Vector<DisjointBoxLayout>&
Realm::getGrids() const
{
  return m_grids;
}

const Vector<ProblemDomain>&
Realm::getDomains() const
{
  return m_domains;
}

Vector<RefCountedPtr<MFLevelGrid>>&
Realm::getMFLevelGrid()
{
  return m_mflg;
}

Vector<RefCountedPtr<MFLevelGrid>>&
Realm::getMFLevelGridCoFi()
{
  return m_mflgCoFi;
}

Vector<RefCountedPtr<MFLevelGrid>>&
Realm::getMFLevelGridFiCo()
{
  return m_mflgFiCo;
}

const RefCountedPtr<EBIndexSpace>&
Realm::getEBIndexSpace(const phase::which_phase a_phase) const
{
  return m_realms[a_phase]->getEBIndexSpace();
}

const Vector<EBISLayout>&
Realm::getEBISLayout(const phase::which_phase a_phase) const
{
  return m_realms[a_phase]->getEBISLayout();
}

const Vector<RefCountedPtr<EBLevelGrid>>&
Realm::getEBLevelGrid(const phase::which_phase a_phase) const
{
  return m_realms[a_phase]->getEBLevelGrid();
}

const Vector<RefCountedPtr<EBLevelGrid>>&
Realm::getEBLevelGridCoFi(const phase::which_phase a_phase) const
{
  return m_realms[a_phase]->getEBLevelGridCoFi();
}

Vector<RefCountedPtr<LayoutData<VoFIterator>>>&
Realm::getVofIterator(const phase::which_phase a_phase) const
{
  return m_realms[a_phase]->getVofIterator();
}

Vector<RefCountedPtr<LayoutData<VoFIterator>>>&
Realm::getMultiCutVofIterator(const phase::which_phase a_phase) const
{
  return m_realms[a_phase]->getMultiCutVofIterator();
}

Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>&
Realm::getFaceIterator(const phase::which_phase a_phase) const
{
  return m_realms[a_phase]->getFaceIterator();
}

Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>&
Realm::getFaceIteratorNoBoundary(const phase::which_phase a_phase) const
{
  return m_realms[a_phase]->getFaceIteratorNoBoundary();
}

Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>&
Realm::getFaceIteratorWithTangentialGhosts(const phase::which_phase a_phase) const
{
  return m_realms[a_phase]->getFaceIteratorWithTangentialGhosts();
}

Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>&
Realm::getMultiCutFaceIterator(const phase::which_phase a_phase) const
{
  return m_realms[a_phase]->getMultiCutFaceIterator();
}

const Vector<RefCountedPtr<EBNonConservativeDivergence>>&
Realm::getNonConservativeDivergence(const phase::which_phase a_phase) const
{
  return m_realms[a_phase]->getNonConservativeDivergence();
}

const Vector<RefCountedPtr<CellCentroidInterpolation>>&
Realm::getCellCentroidInterpolation(const phase::which_phase a_phase) const
{
  return m_realms[a_phase]->getCellCentroidInterpolation();
}

const Vector<RefCountedPtr<EBCentroidInterpolation>>&
Realm::getEBCentroidInterpolation(const phase::which_phase a_phase) const
{
  return m_realms[a_phase]->getEBCentroidInterpolation();
}

Vector<RefCountedPtr<EBCoarAve>>&
Realm::getCoarseAverage(const phase::which_phase a_phase)
{
  return m_realms[a_phase]->getCoarseAverage();
}

EBAMRParticleMesh&
Realm::getParticleMesh(const phase::which_phase a_phase)
{
  return m_realms[a_phase]->getParticleMesh();
}

EBAMRSurfaceDeposition&
Realm::getSurfaceDeposition(const phase::which_phase a_phase)
{
  return m_realms[a_phase]->getSurfaceDeposition();
}

Vector<RefCountedPtr<EBMultigridInterpolator>>&
Realm::getMultigridInterpolator(const phase::which_phase a_phase)
{
  return m_realms[a_phase]->getMultigridInterpolator();
}

const Vector<RefCountedPtr<EBGradient>>&
Realm::getGradientOp(const phase::which_phase a_phase) const
{
  return m_realms[a_phase]->getGradientOp();
}

Vector<RefCountedPtr<EBGhostCellInterpolator>>&
Realm::getGhostCellInterpolator(const phase::which_phase a_phase)
{
  return m_realms[a_phase]->getGhostCellInterpolator();
}

Vector<RefCountedPtr<EBCoarseToFineInterp>>&
Realm::getFineInterp(const phase::which_phase a_phase)
{
  return m_realms[a_phase]->getFineInterp();
}

Vector<RefCountedPtr<EBReflux>>&
Realm::getFluxRegister(const phase::which_phase a_phase)
{
  return m_realms[a_phase]->getFluxRegister();
}

Vector<RefCountedPtr<EBFluxRedistribution>>&
Realm::getRedistributionOp(const phase::which_phase a_phase)
{
  return m_realms[a_phase]->getRedistributionOp();
}

const EBAMRFAB&
Realm::getLevelset(const phase::which_phase a_phase) const
{
  return m_realms[a_phase]->getLevelset();
}

const EBAMRCellData&
Realm::getRegularCells(const phase::which_phase a_phase) const
{
  return m_realms[a_phase]->getRegularCells();
}

const EBAMRCellData&
Realm::getCoveredCells(const phase::which_phase a_phase) const
{
  return m_realms[a_phase]->getCoveredCells();
}

const EBAMRCellData&
Realm::getNotCoveredCells(const phase::which_phase a_phase) const
{
  return m_realms[a_phase]->getNotCoveredCells();
}

const EBAMRCellData&
Realm::getIrregularCells(const phase::which_phase a_phase) const
{
  return m_realms[a_phase]->getIrregularCells();
}

const AMRMask&
Realm::getMask(const std::string& a_mask, const int a_buffer) const
{
  if (!this->queryMask(a_mask, a_buffer)) {
    std::string str = "Realm::getMask - could not find mask '" + a_mask + "'";
    MayDay::Abort(str.c_str());
  }

  return m_masks.at(std::pair<std::string, int>(a_mask, a_buffer));
}

const AMRMask&
Realm::getValidCells() const
{
  return m_validCells;
}

const Vector<RefCountedPtr<LevelTiles>>&
Realm::getLevelTiles() const noexcept
{
  return m_levelTiles;
}

const Vector<RefCountedPtr<MergeBlockGrid>>&
Realm::getMergeBlockGrid() const noexcept
{
  return m_mergeBlockGrids;
}

const AMRParticleGhostMask&
Realm::getParticleGhostMask(const int a_width) const noexcept
{
  const auto it = m_particleGhostMask.find(a_width);
  if (it == m_particleGhostMask.end()) {
    const std::string msg = "Realm::getParticleGhostMask -- width = " + std::to_string(a_width) +
                            " was not registered (call registerParticleGhostMask before regridding)";
    MayDay::Abort(msg.c_str());
  }

  return it->second;
}

const AMRParticleGhostMask&
Realm::getParticleGhostMaskFineToCoar(const int a_width) const noexcept
{
  const auto it = m_particleGhostMaskFineToCoar.find(a_width);
  if (it == m_particleGhostMaskFineToCoar.end()) {
    const std::string msg = "Realm::getParticleGhostMaskFineToCoar -- width = " + std::to_string(a_width) +
                            " was not registered (call registerParticleGhostMask before regridding)";
    MayDay::Abort(msg.c_str());
  }

  return it->second;
}

const AMRParticleGhostMask&
Realm::getParticleGhostMaskCoarToFine(const int a_width) const noexcept
{
  const auto it = m_particleGhostMaskCoarToFine.find(a_width);
  if (it == m_particleGhostMaskCoarToFine.end()) {
    const std::string msg = "Realm::getParticleGhostMaskCoarToFine -- width = " + std::to_string(a_width) +
                            " was not registered (call registerParticleGhostMask before regridding)";
    MayDay::Abort(msg.c_str());
  }

  return it->second;
}

const AMRParticleGhostMask&
Realm::getTrivialParticleGhostMask() const noexcept
{
  static const AMRParticleGhostMask trivial;

  return trivial;
}

Realm::LevelAndBox
Realm::getLevelAndBox(const RealVect& a_pos) const noexcept
{
  // Finest level whose min-tile contains the point wins.
  for (int lvl = m_finestLevel; lvl >= 0; lvl--) {
    IntVect tile;
    for (int dir = 0; dir < SpaceDim; dir++) {
      tile[dir] = static_cast<int>(std::floor((a_pos[dir] - m_probLo[dir]) / (m_minBlockSize * m_dx[lvl])));
    }

    const LevelTiles& tiles = *m_levelTiles[lvl];

    const auto& myTiles = tiles.getMyTiles();
    const auto  mit     = myTiles.find(tile);
    if (mit != myTiles.end()) {
      return LevelAndBox{lvl, mit->second, procID(), true};
    }

#ifdef CH_MPI
    const auto& otherTiles = tiles.getOtherTiles();
    const auto  oit        = otherTiles.find(tile);
    if (oit != otherTiles.end()) {
      return LevelAndBox{lvl, oit->second.first, static_cast<int>(oit->second.second), true};
    }
#endif
  }

  return LevelAndBox{-1, 0u, -1, false};
}

#ifdef CH_USE_PETSC
const RefCountedPtr<PetscGrid>&
Realm::getPetscGrid() const noexcept
{
  return m_petscGrid;
}
#endif

#include <CD_NamespaceFooter.H>
