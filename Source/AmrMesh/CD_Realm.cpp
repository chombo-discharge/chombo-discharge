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
  this->defineGhostTargets();
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
Realm::defineGhostTargets() noexcept
{
  CH_TIME("Realm::defineGhostTargets");
  if (m_verbosity > 5) {
    pout() << "Realm::defineGhostTargets" << endl;
  }

  // Ghosted-box half-width (number of ghost cells). A target box T accepts a scattered particle whose
  // (level-mapped) cell lies in grow(T, ghost), clamped to the domain.
  const int ghost = m_numGhost;

  m_ghostTargetsSame.resize(1 + m_finestLevel);
  m_ghostTargetsCoar.resize(1 + m_finestLevel);
  m_ghostTargetsFine.resize(1 + m_finestLevel);

  // Helper: record, into a per-box source-cell map, every cell of a_cells as scattering to a_target.
  auto registerCells = [](GhostTargetMap& a_map, const Box& a_cells, const GhostTarget& a_target) -> void {
    for (BoxIterator bit(a_cells); bit.ok(); ++bit) {
      a_map[bit()].push_back(a_target);
    }
  };

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl    = m_grids[lvl];
    const Box                domain = m_domains[lvl].domainBox();
    const DataIterator&      dit    = dbl.dataIterator();
    const int                nbox   = dit.size();

    m_ghostTargetsSame[lvl] = RefCountedPtr<LayoutData<GhostTargetMap>>(new LayoutData<GhostTargetMap>(dbl));
    m_ghostTargetsCoar[lvl] = RefCountedPtr<LayoutData<GhostTargetMap>>(new LayoutData<GhostTargetMap>(dbl));
    m_ghostTargetsFine[lvl] = RefCountedPtr<LayoutData<GhostTargetMap>>(new LayoutData<GhostTargetMap>(dbl));

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din    = dit[mybox];
      const Box        srcBox = dbl[din];

      // ---- SAME LEVEL: cells of this box that fall in a same-level neighbour's ghosted box ----
      {
        GhostTargetMap&  map = (*m_ghostTargetsSame[lvl])[din];
        NeighborIterator nit(dbl);
        for (nit.begin(din); nit.ok(); ++nit) {
          const GhostTarget target(dbl.index(nit()), dbl.procID(nit()));

          Box grownNeighbor = grow(dbl[nit()], ghost);
          grownNeighbor &= domain;

          const Box overlap = grownNeighbor & srcBox;
          if (!overlap.isEmpty()) {
            registerCells(map, overlap, target);
          }
        }
      }

      // ---- COARSER LEVEL (lvl -> lvl-1): fine source cells whose coarse parent lies in a coarse box's
      //      ghosted box. Key stays the fine (source-level) cell via refine() of the coarse ghosted box.
      if (lvl > 0) {
        GhostTargetMap&          map        = (*m_ghostTargetsCoar[lvl])[din];
        const DisjointBoxLayout& dblCoar    = m_grids[lvl - 1];
        const Box                domainCoar = m_domains[lvl - 1].domainBox();
        const int                refRat     = m_refinementRatios[lvl - 1];

        for (LayoutIterator lit = dblCoar.layoutIterator(); lit.ok(); ++lit) {
          const GhostTarget target(dblCoar.index(lit()), dblCoar.procID(lit()));

          Box grownCoar = grow(dblCoar[lit()], ghost);
          grownCoar &= domainCoar;

          const Box overlap = refine(grownCoar, refRat) & srcBox; // fine cells covered by the coarse ghosted box
          if (!overlap.isEmpty()) {
            registerCells(map, overlap, target);
          }
        }
      }

      // ---- FINER LEVEL (lvl -> lvl+1): coarse source cells whose refinement overlaps a fine box's
      //      ghosted box. Key stays the coarse (source-level) cell via coarsen() of the fine ghosted box.
      if (lvl < m_finestLevel) {
        GhostTargetMap&          map        = (*m_ghostTargetsFine[lvl])[din];
        const DisjointBoxLayout& dblFine    = m_grids[lvl + 1];
        const Box                domainFine = m_domains[lvl + 1].domainBox();
        const int                refRat     = m_refinementRatios[lvl];

        for (LayoutIterator lit = dblFine.layoutIterator(); lit.ok(); ++lit) {
          const GhostTarget target(dblFine.index(lit()), dblFine.procID(lit()));

          Box grownFine = grow(dblFine[lit()], ghost);
          grownFine &= domainFine;

          const Box overlap = coarsen(grownFine, refRat) & srcBox; // coarse cells overlapping the fine ghosted box
          if (!overlap.isEmpty()) {
            registerCells(map, overlap, target);
          }
        }
      }
    }
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

const AMRGhostTargets&
Realm::getGhostTargetsSame() const noexcept
{
  return m_ghostTargetsSame;
}

const AMRGhostTargets&
Realm::getGhostTargetsCoar() const noexcept
{
  return m_ghostTargetsCoar;
}

const AMRGhostTargets&
Realm::getGhostTargetsFine() const noexcept
{
  return m_ghostTargetsFine;
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
