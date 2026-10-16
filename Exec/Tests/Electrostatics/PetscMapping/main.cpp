/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   main.cpp
 * @brief  Tests of the mapping between chombo-discharge and PETSc.
 * @author Robert Marskar
 */

// Std includes
#include <random>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_Driver.H>
#include <CD_FieldSolverGMG.H>
#include <CD_RodDielectric.H>
#include <CD_FieldStepper.H>
#include <CD_ParallelOps.H>

using namespace ChomboDischarge;
using namespace Physics::Electrostatics;

#ifdef CH_USE_PETSC
/**
 * @brief Fill every cell of a data holder, ghost cells included, with uniform random numbers in [0, 1).
 * @param[in,out] a_data Data holder.
 */
static void
fillRandom(MFAMRCellData& a_data)
{
  std::mt19937_64                  rng(1234 + procID());
  std::uniform_real_distribution<> dist(0.0, 1.0);

  for (int lvl = 0; lvl < a_data.size(); lvl++) {
    for (DataIterator dit = a_data[lvl]->dataIterator(); dit.ok(); ++dit) {
      for (int iphase = 0; iphase < (*a_data[lvl])[dit()].numPhases(); iphase++) {
        FArrayBox& fab = (*a_data[lvl])[dit()].getPhase(iphase).getFArrayBox();

        for (BoxIterator bit(fab.box()); bit.ok(); ++bit) {
          fab(bit(), 0) = dist(rng);
        }
      }
    }
  }
}

/**
 * @brief Abort unless the rows of a stencil are unique and ascending, as PetscGrid returns them.
 * @param[in] a_stencil Stencil to check.
 */
static void
checkCompressed(const PetscStencil& a_stencil)
{
  for (int i = 1; i < a_stencil.size(); i++) {
    if (a_stencil.row(i) <= a_stencil.row(i - 1)) {
      MayDay::Abort("PetscMapping -- stencil rows are not unique and ascending");
    }
  }
}

/**
 * @brief Abort unless PetscStencil::compress sorts the rows and sums the weights of duplicate rows.
 */
static void
testCompress()
{
  PetscStencil a;
  a.add(7, 1.0);
  a.add(3, 2.0);
  a.add(7, 4.0);

  PetscStencil b;
  b.add(3, 8.0);
  b.add(5, 16.0);

  a += b;
  a *= 0.5;
  a.compress();

  checkCompressed(a);

  if (a.size() != 3 || a.row(0) != 3 || a.row(1) != 5 || a.row(2) != 7 || a.weight(0) != 5.0 || a.weight(1) != 8.0 ||
      a.weight(2) != 2.5) {
    MayDay::Abort("PetscMapping -- PetscStencil::compress did not merge duplicate rows");
  }
}

/**
 * @brief Run the tests.
 * @details Aborts on the first failure.
 * - PetscStencil::compress sorts the rows and merges duplicates, and every stencil from PetscGrid is compressed.
 * - The Chombo -> PETSc -> Chombo round trip reproduces every row exactly.
 * - Ghost cells that overlap a neighbouring box hold the global row of that box's cell.
 * - PetscGrid::getPetscStencilWithFinerLevel finds the level l+1 rows under a covered level l cell: their average
 *   equals the conservative average.
 * - Every coarse-fine ghost cell that the multigrid interpolator fills: its stencil from
 *   EBMultigridInterpolator::getGhostStencil, mapped by PetscGrid::getPetscStencilWithCoarserLevel and applied to the
 *   PETSc vector, reproduces the value from AmrMesh::interpGhostMG to round-off.
 * @param[in] a_amr   AmrMesh, after the driver has set up the grids.
 * @param[in] a_realm Realm.
 */
static void
runTests(const AmrMesh& a_amr, const std::string& a_realm)
{
  const RefCountedPtr<PetscGrid>& petscGrid   = a_amr.getPetscGrid(a_realm);
  const int                       finestLevel = a_amr.getFinestLevel();
  const int                       numPhases   = petscGrid->getNumPhases();

  const Vector<RefCountedPtr<LevelData<BaseFab<PetscAMRCell>>>>& amrToPetsc     = petscGrid->getAMRToPetsc();
  const Vector<RefCountedPtr<MFLevelGrid>>&                      mflg           = petscGrid->getMFLevelGrids();
  const Vector<RefCountedPtr<LevelData<BaseFab<PetscAMRCell>>>>& amrToPetscFiCo = petscGrid->getAMRToPetscFiCo();

  if (finestLevel < 1) {
    MayDay::Abort("PetscMapping -- the test needs at least two levels");
  }

  MFAMRCellData phi;
  MFAMRCellData phiBack;

  a_amr.allocate(phi, a_realm, 1);
  a_amr.allocate(phiBack, a_realm, 1);

  fillRandom(phi);

  // Make the data on each phase consistent with a composite solution: covered cells hold the conservative average
  // and ghost cells overlapping neighbouring boxes hold their values.
  for (int iphase = 0; iphase < numPhases; iphase++) {
    EBAMRCellData phiPhase = a_amr.alias(phase::which_phase(iphase), phi);

    a_amr.conservativeAverage(phiPhase, a_realm, phase::which_phase(iphase));

    phiPhase.exchange();
  }

  Vec x;
  Vec xAll;

  VecScatter toAll;

  petscGrid->create(x);
  petscGrid->putChomboInPetsc(x, phi);

  // Every rank gets the whole vector, so that ghost cells can be checked against rows owned by other ranks.
  PetscCallVoid(VecScatterCreateToAll(x, &toAll, &xAll));
  PetscCallVoid(VecScatterBegin(toAll, x, xAll, INSERT_VALUES, SCATTER_FORWARD));
  PetscCallVoid(VecScatterEnd(toAll, x, xAll, INSERT_VALUES, SCATTER_FORWARD));

  const PetscScalar* xAllArr;
  PetscCallVoid(VecGetArrayRead(xAll, &xAllArr));

  testCompress();

  // 1. Round trip.
  DataOps::setValue(phiBack, 0.0);

  petscGrid->putPetscInChombo(phiBack, x);

  long long numRoundTrip = 0;

  for (int lvl = 0; lvl <= finestLevel; lvl++) {
    for (DataIterator dit = phi[lvl]->dataIterator(); dit.ok(); ++dit) {
      const BaseFab<PetscAMRCell>& cells = (*amrToPetsc[lvl])[dit()];

      for (int iphase = 0; iphase < numPhases; iphase++) {
        const FArrayBox& orig = (*phi[lvl])[dit()].getPhase(iphase).getFArrayBox();
        const FArrayBox& back = (*phiBack[lvl])[dit()].getPhase(iphase).getFArrayBox();

        for (BoxIterator bit(mflg[lvl]->getGrids()[dit()]); bit.ok(); ++bit) {
          if (cells(bit()).getPetscRow(iphase) >= 0) {
            if (orig(bit(), 0) != back(bit(), 0)) {
              MayDay::Abort("PetscMapping -- round trip does not reproduce the data");
            }

            numRoundTrip++;
          }
        }
      }
    }
  }

  // 2. Ghost cells overlapping neighbouring boxes hold the neighbour's global row.
  long long numNeighborGhosts = 0;

  for (int lvl = 0; lvl <= finestLevel; lvl++) {
    for (DataIterator dit = phi[lvl]->dataIterator(); dit.ok(); ++dit) {
      const BaseFab<PetscAMRCell>& cells = (*amrToPetsc[lvl])[dit()];
      const Box                    box   = mflg[lvl]->getGrids()[dit()];

      for (int iphase = 0; iphase < numPhases; iphase++) {
        const FArrayBox& data = (*phi[lvl])[dit()].getPhase(iphase).getFArrayBox();

        for (BoxIterator bit(cells.box()); bit.ok(); ++bit) {
          const PetscInt row = cells(bit()).getPetscRow(iphase);

          if (!box.contains(bit()) && row >= 0) {
            if (data(bit(), 0) != xAllArr[row]) {
              MayDay::Abort("PetscMapping -- ghost cell does not hold the global row of the neighbouring cell");
            }

            numNeighborGhosts++;
          }
        }
      }
    }
  }

  // 3. Two-level stencils: the average of the level l+1 rows under a regular covered cell on level l equals the
  // conservative average, where those level l+1 cells are rows.
  long long numCovered = 0;

  for (int lvl = 0; lvl < finestLevel; lvl++) {
    const int refRat = a_amr.getRefinementRatios()[lvl];

    for (DataIterator dit = phi[lvl]->dataIterator(); dit.ok(); ++dit) {
      const BaseFab<PetscAMRCell>& cells = (*amrToPetsc[lvl])[dit()];
      const Box                    box   = mflg[lvl]->getGrids()[dit()];

      for (int iphase = 0; iphase < numPhases; iphase++) {
        const EBISBox&   ebisBox = mflg[lvl]->getEBLevelGrid(iphase).getEBISL()[dit()];
        const FArrayBox& data    = (*phi[lvl])[dit()].getPhase(iphase).getFArrayBox();

        for (BoxIterator bit(box); bit.ok(); ++bit) {
          if (cells(bit()).isCoveredByFinerGrid() && ebisBox.isRegular(bit())) {
            const Box children = refine(Box(bit(), bit()), refRat);

            // Skip cells whose children are themselves covered by level l+2.
            const BaseFab<PetscAMRCell>& fineCells = (*amrToPetscFiCo[lvl + 1])[dit()];

            bool childrenValid = true;
            for (BoxIterator cit(children); cit.ok(); ++cit) {
              childrenValid = childrenValid && !fineCells(cit()).isCoveredByFinerGrid();
            }
            if (!childrenValid) {
              continue;
            }

            VoFStencil fineStencil;

            for (BoxIterator cit(children); cit.ok(); ++cit) {
              fineStencil.add(VolIndex(cit(), 0), 1.0 / children.numPts());
            }

            const PetscStencil stencil = petscGrid->getPetscStencilWithFinerLevel(VoFStencil(),
                                                                                  fineStencil,
                                                                                  lvl,
                                                                                  dit(),
                                                                                  iphase);

            checkCompressed(stencil);

            Real value = 0.0;
            for (int i = 0; i < stencil.size(); i++) {
              value += stencil.weight(i) * xAllArr[stencil.row(i)];
            }

            if (std::abs(value - data(bit(), 0)) > 1.E-12) {
              MayDay::Abort("PetscMapping -- average of fine rows does not equal the conservative average");
            }

            numCovered++;
          }
        }
      }
    }
  }

  // 4. Coarse-fine ghost cells.
  long long numCFGhosts = 0;
  Real      maxCFError  = 0.0;

  for (int iphase = 0; iphase < numPhases; iphase++) {
    const phase::which_phase curPhase = phase::which_phase(iphase);

    EBAMRCellData phiPhase = a_amr.alias(curPhase, phi);

    a_amr.interpGhostMG(phiPhase, a_realm, curPhase);

    const Vector<RefCountedPtr<EBMultigridInterpolator>>& interpolators = a_amr.getMultigridInterpolator(a_realm,
                                                                                                         curPhase);

    for (int lvl = 1; lvl <= finestLevel; lvl++) {
      const EBLevelGrid&       eblg    = mflg[lvl]->getEBLevelGrid(iphase);
      const DisjointBoxLayout& dbl     = eblg.getDBL();
      const ProblemDomain&     domain  = eblg.getDomain();
      const int                ghostCF = interpolators[lvl]->getGhostCF();

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
        const Box                    box     = dbl[dit()];
        const EBISBox&               ebisBox = eblg.getEBISL()[dit()];
        const BaseFab<PetscAMRCell>& cells   = (*amrToPetsc[lvl])[dit()];
        const EBCellFAB&             data    = (*phiPhase[lvl])[dit()];

        // The ghost cells coarseFineInterp fills: the first layer across each face, and cells within ghostCF of the
        // cut cells in this box.
        IntVectSet candidates;

        for (int dir = 0; dir < SpaceDim; dir++) {
          for (SideIterator sit; sit.ok(); ++sit) {
            candidates |= adjCellBox(box, dir, sit(), 1);
          }
        }

        IntVectSet nearEB = ebisBox.getIrregIVS(box);
        nearEB.grow(ghostCF);
        candidates |= nearEB;

        for (IVSIterator ivsit(candidates); ivsit.ok(); ++ivsit) {
          const IntVect iv = ivsit();

          if (box.contains(iv) || !domain.contains(iv) || !cells(iv).isGhostCF() || ebisBox.isCovered(iv)) {
            continue;
          }

          const VolIndex vof(iv, 0);

          // The ghost cell's value is what the interpolator fills it with: a stencil over level lvl cells and the
          // coarsened level lvl-1 cells under this grid.
          const std::pair<VoFStencil, VoFStencil> ghost = interpolators[lvl]->getGhostStencil(vof, dit());

          const PetscStencil stencil = petscGrid->getPetscStencilWithCoarserLevel(ghost.first,
                                                                                  ghost.second,
                                                                                  lvl,
                                                                                  dit(),
                                                                                  iphase);

          checkCompressed(stencil);

          Real value = 0.0;
          Real scale = 0.0;

          for (int i = 0; i < stencil.size(); i++) {
            value += stencil.weight(i) * xAllArr[stencil.row(i)];
            scale += std::abs(stencil.weight(i));
          }

          const Real error = std::abs(value - data(vof, 0));

          maxCFError = std::max(maxCFError, error);

          if (error > 1.E-10 * std::max(1.0, scale)) {
            pout() << "PetscMapping -- level = " << lvl << ", phase = " << iphase << ", ghost cell = " << iv
                   << ", interpGhostMG = " << data(vof, 0) << ", PETSc stencil = " << value << endl;

            MayDay::Abort("PetscMapping -- PETSc ghost cell stencil does not reproduce interpGhostMG");
          }

          numCFGhosts++;
        }
      }
    }
  }

  PetscCallVoid(VecRestoreArrayRead(xAll, &xAllArr));
  PetscCallVoid(VecScatterDestroy(&toAll));
  PetscCallVoid(VecDestroy(&xAll));

  petscGrid->destroy(x);

  numRoundTrip      = ParallelOps::sum(numRoundTrip);
  numNeighborGhosts = ParallelOps::sum(numNeighborGhosts);
  numCovered        = ParallelOps::sum(numCovered);
  numCFGhosts       = ParallelOps::sum(numCFGhosts);
  maxCFError        = ParallelOps::max(maxCFError);

  pout() << "PetscMapping -- round trip: " << numRoundTrip << " rows" << endl;
  pout() << "PetscMapping -- neighbour ghost cells: " << numNeighborGhosts << endl;
  pout() << "PetscMapping -- covered cells: " << numCovered << endl;
  pout() << "PetscMapping -- coarse-fine ghost cells: " << numCFGhosts << ", max error = " << maxCFError << endl;

  if (numRoundTrip == 0 || numNeighborGhosts == 0 || numCovered == 0 || numCFGhosts == 0) {
    MayDay::Abort("PetscMapping -- a test checked no cells");
  }

  pout() << "PetscMapping -- all tests passed" << endl;
}
#endif

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

#ifdef CH_USE_PETSC
  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new RodDielectric());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto tagger      = RefCountedPtr<CellTagger>(nullptr);
  auto timestepper = RefCountedPtr<FieldStepper<FieldSolverGMG>>(new FieldStepper<FieldSolverGMG>());
  auto engine      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

  engine->setupAndRun();

  std::string realm;

  ParmParse pp("FieldStepper");
  pp.get("realm", realm);

  runTests(*amr, realm);
#else
  MayDay::Abort("PetscMapping -- this test requires USE_PETSC=TRUE");
#endif

  ChomboDischarge::finalize();
}
