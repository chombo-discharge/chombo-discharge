/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   main.cpp
 * @brief  Electrostatics in the PartialDischarge geometry (a gas void in a dielectric), with error diagnostics.
 * @author Robert Marskar
 */

// Std includes
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_Driver.H>
#include <CD_Aerosol.H>
#include <CD_RefinementBoxTagger.H>
#include <CD_FieldStepper.H>
#include <CD_FieldSolverGMG.H>
#ifdef CH_USE_PETSC
#include <CD_FieldSolverAMG.H>
#endif
#include <CD_Location.H>
#include <CD_ParallelOps.H>

using namespace ChomboDischarge;
using namespace Physics::Electrostatics;

/**
 * @brief Exact potential of a gas void (sphere in 3D, cylinder in 2D) in a dielectric, in the uniform field set up by
 * Dirichlet conditions phi = 0 and phi = 1 on the low and high y sides.
 * @details Inside the void the field is uniform. Outside it is the uniform field plus a dipole. The void is small
 * compared to the domain, so the images of the dipole in the domain sides are neglected.
 */
class ExactSolution
{
public:
  /**
   * @brief Constructor.
   * @param[in] a_center      Void center
   * @param[in] a_radius      Void radius
   * @param[in] a_epsGas      Permittivity of the gas (inside the void)
   * @param[in] a_epsSolid    Permittivity of the dielectric (outside the void)
   * @param[in] a_loY         Low y coordinate of the domain, where phi = 0
   * @param[in] a_hiY         High y coordinate of the domain, where phi = 1
   */
  ExactSolution(const RealVect& a_center,
                const Real      a_radius,
                const Real      a_epsGas,
                const Real      a_epsSolid,
                const Real      a_loY,
                const Real      a_hiY)
    : m_center(a_center),
      m_radius(a_radius),
      m_grad(1.0 / (a_hiY - a_loY)),
      m_phiC((a_center[1] - a_loY) / (a_hiY - a_loY))
  {

    const Real ei = a_epsGas;
    const Real eo = a_epsSolid;

    if (SpaceDim == 3) {
      m_inside  = 3.0 * eo / (ei + 2.0 * eo);
      m_outside = (ei - eo) / (ei + 2.0 * eo);
    }
    else {
      m_inside  = 2.0 * eo / (ei + eo);
      m_outside = (ei - eo) / (ei + eo);
    }
  }

  /**
   * @brief Exact potential.
   * @param[in] a_pos   Position
   * @param[in] a_phase Phase whose solution is evaluated (extended analytically outside the phase)
   * @return The potential
   */
  Real
  phi(const RealVect& a_pos, const int a_phase) const
  {
    const RealVect x = a_pos - m_center;
    const Real     y = x[1];

    if (a_phase == phase::gas) {
      return m_phiC + m_grad * m_inside * y;
    }

    const Real r = std::max(x.vectorLength(), 1.E-30);

    return m_phiC + m_grad * y * (1.0 - m_outside * std::pow(m_radius / r, SpaceDim));
  }

  /**
   * @brief Exact electric field, E = -grad(phi).
   * @param[in] a_pos   Position
   * @param[in] a_phase Phase whose solution is evaluated (extended analytically outside the phase)
   * @return The electric field
   */
  RealVect
  E(const RealVect& a_pos, const int a_phase) const
  {
    RealVect grad = RealVect::Zero;

    const RealVect x = a_pos - m_center;
    const Real     y = x[1];

    if (a_phase == phase::gas) {
      grad[1] = m_grad * m_inside;
    }
    else {
      const Real r  = std::max(x.vectorLength(), 1.E-30);
      const Real g  = 1.0 - m_outside * std::pow(m_radius / r, SpaceDim);
      const Real dg = SpaceDim * m_outside * std::pow(m_radius, SpaceDim) / std::pow(r, SpaceDim + 1);

      for (int dir = 0; dir < SpaceDim; dir++) {
        grad[dir] = m_grad * ((dir == 1 ? g : 0.0) + y * dg * x[dir] / r);
      }
    }

    return -grad;
  }

  /**
   * @brief Background field strength |E_0|.
   * @return |E_0|
   */
  Real
  E0() const
  {
    return m_grad;
  }

  /**
   * @brief Signed distance to the void surface (positive outside the void).
   * @param[in] a_pos Position
   * @return The distance
   */
  Real
  distance(const RealVect& a_pos) const
  {
    return (a_pos - m_center).vectorLength() - m_radius;
  }

  /**
   * @brief Void radius.
   * @return The radius
   */
  Real
  radius() const
  {
    return m_radius;
  }

protected:
  /**
   * @brief Void center.
   */
  RealVect m_center;

  /**
   * @brief Void radius.
   */
  Real m_radius;

  /**
   * @brief Gradient of the background potential, 1/L_y.
   */
  Real m_grad;

  /**
   * @brief Background potential at the void center.
   */
  Real m_phiC;

  /**
   * @brief Field enhancement inside the void.
   */
  Real m_inside;

  /**
   * @brief Dipole coefficient outside the void.
   */
  Real m_outside;
};

/**
 * @brief FieldStepper that compares the solution and its gradient with the exact solution.
 */
template <class T>
class DielectricVoidStepper : public FieldStepper<T>
{
public:
  /**
   * @brief Compare the solution, its gradient and the centroid-interpolated gradient with the exact solution.
   * @details Reports, per phase, level and cell type, the errors of
   *  A: the computed potential at cell centers;
   *  B: E from the computed potential, at cell centers;
   *  C: E from the exact potential at cell centers (the gradient operator alone), at cell centers;
   *  D: B interpolated to cell centroids (what the plot files show), against the exact E at the centroids;
   *  E: the exact E at cell centers interpolated to the centroids (the interpolation alone);
   *  F: C interpolated to the centroids (gradient and interpolation, no solver).
   * Only valid cells (not covered by a finer level) within a_maxRadius of the void center are included. The E errors
   * are relative to |E_0|, the potential errors to |E_0| times the void radius.
   * @param[in] a_exact     Exact solution
   * @param[in] a_maxRadius Largest distance from the void center of the cells that are included
   * @param[in] a_numWorst  Number of worst cut cells to list, per phase
   */
  void
  diagnose(const ExactSolution& a_exact, const Real a_maxRadius, const int a_numWorst) const
  {
    const RefCountedPtr<AmrMesh>&     amr    = this->m_amr;
    const RefCountedPtr<FieldSolver>& solver = this->m_fieldSolver;
    const std::string&                realm  = this->m_realm;

    const int           finestLevel = amr->getFinestLevel();
    const RealVect      probLo      = amr->getProbLo();
    const Vector<Real>& dx          = amr->getDx();

    // Exact potential and exact E at the cell centers, in all cells including ghost cells.
    MFAMRCellData phiExact;
    MFAMRCellData exactCenter;

    amr->allocate(phiExact, realm, 1);
    amr->allocate(exactCenter, realm, SpaceDim);

    for (int lvl = 0; lvl <= finestLevel; lvl++) {
      const DataIterator& dit  = amr->getGrids(realm)[lvl].dataIterator();
      const int           nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        for (int iphase = 0; iphase < (*phiExact[lvl])[din].numPhases(); iphase++) {
          EBCellFAB&     phi     = (*phiExact[lvl])[din].getPhase(iphase);
          EBCellFAB&     E       = (*exactCenter[lvl])[din].getPhase(iphase);
          const EBISBox& ebisbox = phi.getEBISBox();

          phi.setVal(0.0);
          E.setVal(0.0);

          for (VoFIterator vofit(IntVectSet(phi.box()), ebisbox.getEBGraph()); vofit.ok(); ++vofit) {
            const VolIndex& vof = vofit();
            const RealVect  pos = probLo + Location::position(Location::Cell::Center, vof, ebisbox, dx[lvl]);
            const RealVect  e   = a_exact.E(pos, iphase);

            phi(vof, 0) = a_exact.phi(pos, iphase);

            for (int dir = 0; dir < SpaceDim; dir++) {
              E(vof, dir) = e[dir];
            }
          }
        }
      }
    }

    // B: E from the computed potential. C: E from the exact potential.
    MFAMRCellData computedCenter;
    MFAMRCellData gradExactCenter;

    amr->allocate(computedCenter, realm, SpaceDim);
    amr->allocate(gradExactCenter, realm, SpaceDim);

    solver->computeElectricField(computedCenter, solver->getPotential());
    solver->computeElectricField(gradExactCenter, phiExact);

    // D, E, F: the same fields interpolated to the cell centroids, as when the plot files are written.
    MFAMRCellData computedCentroid;
    MFAMRCellData exactInterpCentroid;
    MFAMRCellData gradExactCentroid;

    amr->allocate(computedCentroid, realm, SpaceDim);
    amr->allocate(exactInterpCentroid, realm, SpaceDim);
    amr->allocate(gradExactCentroid, realm, SpaceDim);

    amr->copyData(computedCentroid, computedCenter);
    amr->copyData(exactInterpCentroid, exactCenter);
    amr->copyData(gradExactCentroid, gradExactCenter);

    // A rank can own no boxes on a level, so the number of phases is taken over all levels and ranks.
    int numPhases = 0;

    for (int lvl = 0; lvl <= finestLevel; lvl++) {
      const DataIterator& dit  = amr->getGrids(realm)[lvl].dataIterator();
      const int           nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(max : numPhases)
      for (int mybox = 0; mybox < nbox; mybox++) {
        numPhases = std::max(numPhases, (*phiExact[lvl])[dit[mybox]].numPhases());
      }
    }

    numPhases = ParallelOps::max(numPhases);

    for (int iphase = 0; iphase < numPhases; iphase++) {
      const phase::which_phase p = (iphase == 0) ? phase::gas : phase::solid;

      for (MFAMRCellData* data : {&computedCentroid, &exactInterpCentroid, &gradExactCentroid}) {
        EBAMRCellData alias = amr->alias(p, *data);

        amr->interpToCentroids(alias, realm, p);
      }
    }

    // Accumulate the errors.
    constexpr int numMeasures        = 6;
    const char*   names[numMeasures] = {"A phi (solver)",
                                        "B E center",
                                        "C E(exact phi) center",
                                        "D E centroid (plot)",
                                        "E interp(exact E)",
                                        "F interp(C)"};

    const Real phiScale = a_exact.E0() * a_exact.radius();
    const Real eScale   = a_exact.E0();

    const AMRMask& validCells = amr->getValidCells(realm);

    // One record per cut cell: measure D first, so that the records sort by the error seen in the plots.
    constexpr int     recordSize = 7 + SpaceDim + numMeasures;
    std::vector<Real> records;

    pout() << endl
           << "DielectricVoid -- errors near the void (r < " << a_maxRadius << "), relative to |E_0| = " << eScale
           << " (E) and |E_0|*radius (phi)" << endl;

    for (int iphase = 0; iphase < numPhases; iphase++) {
      for (int lvl = 0; lvl <= finestLevel; lvl++) {
        // Index t*numMeasures + m, with t = 0 for regular and t = 1 for cut cells.
        Real      maxErr[2 * numMeasures];
        Real      sumSq[2 * numMeasures];
        long long count[2] = {0, 0};

        for (int i = 0; i < 2 * numMeasures; i++) {
          maxErr[i] = 0.0;
          sumSq[i]  = 0.0;
        }

        const DataIterator& dit  = amr->getGrids(realm)[lvl].dataIterator();
        const int           nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(max : maxErr[ : 2 * numMeasures]) \
  reduction(+ : sumSq[ : 2 * numMeasures], count[ : 2])
        for (int mybox = 0; mybox < nbox; mybox++) {
          const DataIndex&     din     = dit[mybox];
          const Box            box     = amr->getGrids(realm)[lvl][din];
          const BaseFab<bool>& valid   = (*validCells[lvl])[din];
          const EBCellFAB&     phiH    = (*solver->getPotential()[lvl])[din].getPhase(iphase);
          const EBCellFAB&     phiEx   = (*phiExact[lvl])[din].getPhase(iphase);
          const EBCellFAB&     Bfab    = (*computedCenter[lvl])[din].getPhase(iphase);
          const EBCellFAB&     Cfab    = (*gradExactCenter[lvl])[din].getPhase(iphase);
          const EBCellFAB&     Dfab    = (*computedCentroid[lvl])[din].getPhase(iphase);
          const EBCellFAB&     Efab    = (*exactInterpCentroid[lvl])[din].getPhase(iphase);
          const EBCellFAB&     Ffab    = (*gradExactCentroid[lvl])[din].getPhase(iphase);
          const EBCellFAB&     exactC  = (*exactCenter[lvl])[din].getPhase(iphase);
          const EBISBox&       ebisbox = phiH.getEBISBox();

          std::vector<Real> boxRecords;

          for (VoFIterator vofit(IntVectSet(box), ebisbox.getEBGraph()); vofit.ok(); ++vofit) {
            const VolIndex& vof = vofit();
            const IntVect&  iv  = vof.gridIndex();

            if (!valid(iv)) {
              continue;
            }

            const RealVect center   = probLo + Location::position(Location::Cell::Center, vof, ebisbox, dx[lvl]);
            const RealVect centroid = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx[lvl]);

            if (a_exact.distance(center) + a_exact.radius() > a_maxRadius) {
              continue;
            }

            const RealVect exactE    = RealVect(D_DECL(exactC(vof, 0), exactC(vof, 1), exactC(vof, 2)));
            const RealVect exactEcen = a_exact.E(centroid, iphase);

            auto diff = [&](const EBCellFAB& a_fab, const RealVect& a_ref) -> Real {
              RealVect d;

              for (int dir = 0; dir < SpaceDim; dir++) {
                d[dir] = a_fab(vof, dir) - a_ref[dir];
              }

              return d.vectorLength() / eScale;
            };

            Real err[numMeasures];

            err[0] = std::abs(phiH(vof, 0) - phiEx(vof, 0)) / phiScale;
            err[1] = diff(Bfab, exactE);
            err[2] = diff(Cfab, exactE);
            err[3] = diff(Dfab, exactEcen);
            err[4] = diff(Efab, exactEcen);
            err[5] = diff(Ffab, exactEcen);

            const int t = ebisbox.isIrregular(iv) ? 1 : 0;

            count[t]++;

            for (int m = 0; m < numMeasures; m++) {
              maxErr[t * numMeasures + m] = std::max(maxErr[t * numMeasures + m], err[m]);
              sumSq[t * numMeasures + m] += err[m] * err[m];
            }

            if (t == 1) {
              boxRecords.push_back(err[3]);
              boxRecords.push_back(iphase);
              boxRecords.push_back(lvl);
              boxRecords.push_back(ebisbox.volFrac(vof));
              boxRecords.push_back(a_exact.distance(centroid) / dx[lvl]);
              boxRecords.push_back(ebisbox.bndryArea(vof));
              boxRecords.push_back(ebisbox.numVoFs(iv));

              for (int dir = 0; dir < SpaceDim; dir++) {
                boxRecords.push_back(iv[dir]);
              }

              for (int m = 0; m < numMeasures; m++) {
                boxRecords.push_back(err[m]);
              }
            }
          }

#pragma omp critical(DielectricVoidRecords)
          records.insert(records.end(), boxRecords.begin(), boxRecords.end());
        }

        for (int t = 0; t < 2; t++) {
          const long long n = ParallelOps::sum(count[t]);

          if (n == 0) {
            continue;
          }

          pout() << "  phase " << iphase << " level " << lvl << (t == 1 ? " cut    " : " regular") << " cells " << n
                 << endl;

          for (int m = 0; m < numMeasures; m++) {
            const Real mx  = ParallelOps::max(maxErr[t * numMeasures + m]);
            const Real rms = std::sqrt(ParallelOps::sum(sumSq[t * numMeasures + m]) / n);

            pout() << "      " << names[m] << ": max " << mx << ", rms " << rms << endl;
          }
        }
      }
    }

    // The worst cut cells by measure D, gathered over ranks.
    std::vector<Real> allRecords = records;

#ifdef CH_MPI
    int localSize = records.size();

    std::vector<int> sizes(numProc());

    MPI_Allgather(&localSize, 1, MPI_INT, sizes.data(), 1, MPI_INT, Chombo_MPI::comm);

    std::vector<int> displs(numProc(), 0);

    for (int i = 1; i < numProc(); i++) {
      displs[i] = displs[i - 1] + sizes[i - 1];
    }

    allRecords.resize(displs.back() + sizes.back());

    MPI_Allgatherv(records.data(),
                   localSize,
                   MPI_CH_REAL,
                   allRecords.data(),
                   sizes.data(),
                   displs.data(),
                   MPI_CH_REAL,
                   Chombo_MPI::comm);
#endif

    const int numRecords = allRecords.size() / recordSize;

    std::vector<int> order(numRecords);

    for (int i = 0; i < numRecords; i++) {
      order[i] = i;
    }

    std::sort(order.begin(), order.end(), [&](const int a, const int b) -> bool {
      return allRecords[a * recordSize] > allRecords[b * recordSize];
    });

    for (int iphase = 0; iphase < numPhases; iphase++) {
      pout() << "  worst cut cells by D, phase " << iphase << ":" << endl;

      int listed = 0;

      for (int i = 0; i < numRecords && listed < a_numWorst; i++) {
        const Real* r = &allRecords[order[i] * recordSize];

        if (static_cast<int>(r[1]) != iphase) {
          continue;
        }

        pout() << "    lvl " << static_cast<int>(r[2]) << " iv (";

        for (int dir = 0; dir < SpaceDim; dir++) {
          pout() << static_cast<int>(r[7 + dir]) << (dir < SpaceDim - 1 ? "," : ")");
        }

        pout() << " kappa " << r[3] << " d/dx " << r[4] << " area " << r[5] << " nvofs " << static_cast<int>(r[6])
               << " |";

        for (int m = 0; m < numMeasures; m++) {
          pout() << " " << names[m][0] << "=" << r[7 + SpaceDim + m];
        }

        pout() << endl;

        listed++;
      }
    }
  }
};

/**
 * @brief Run the solver and the diagnostics.
 * @param[in] a_exact     Exact solution
 * @param[in] a_maxRadius Largest distance from the void center of the cells in the diagnostics
 * @param[in] a_numWorst  Number of worst cells to list
 * @param[in] a_boxTagger Refine with RefinementBoxTagger in addition to the geometric refinement, or not
 */
template <class T>
void
run(const ExactSolution& a_exact, const Real a_maxRadius, const int a_numWorst, const bool a_boxTagger)
{
  std::string realm;

  ParmParse ppStepper("FieldStepper");

  ppStepper.get("realm", realm);

  // The time stepper is owned by 'timestepper'; 'stepper' only refers to it, for the diagnostics.
  DielectricVoidStepper<T>* stepper = new DielectricVoidStepper<T>();

  auto timestepper = RefCountedPtr<TimeStepper>(stepper);
  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new Aerosol());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto tagger      = a_boxTagger ? RefCountedPtr<CellTagger>(new RefinementBoxTagger(amr, realm))
                                 : RefCountedPtr<CellTagger>(nullptr);
  auto driver      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

  driver->setupAndRun();

  stepper->diagnose(a_exact, a_maxRadius, a_numWorst);
}

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  {
    std::string fieldSolver;
    int         numWorst;
    Real        maxRadiusFactor;
    bool        boxTagger;

    ParmParse ppVoid("DielectricVoid");

    ppVoid.get("field_solver", fieldSolver);
    ppVoid.get("num_worst", numWorst);
    ppVoid.get("max_radius", maxRadiusFactor);
    ppVoid.get("box_tagger", boxTagger);

    // The exact solution, from the geometry and domain inputs. The void must be a hole (Aerosol.invert = true).
    ParmParse ppAerosol("Aerosol");
    ParmParse ppSphere("Aerosol.sphere1");
    ParmParse ppAmr("AmrMesh");

    bool         invert;
    Real         epsGas;
    Real         epsSolid;
    Real         radius;
    Vector<Real> center(SpaceDim);
    Vector<Real> lo(SpaceDim);
    Vector<Real> hi(SpaceDim);

    ppAerosol.get("invert", invert);
    ppAerosol.get("eps0", epsGas);
    ppAerosol.get("permittivity", epsSolid);
    ppSphere.get("radius", radius);
    ppSphere.getarr("center", center, 0, SpaceDim);
    ppAmr.getarr("lo_corner", lo, 0, SpaceDim);
    ppAmr.getarr("hi_corner", hi, 0, SpaceDim);

    if (!invert) {
      MayDay::Abort("DielectricVoid -- the exact solution needs Aerosol.invert = true (a gas void)");
    }

    const ExactSolution exact(RealVect(D_DECL(center[0], center[1], center[2])),
                              radius,
                              epsGas,
                              epsSolid,
                              lo[1],
                              hi[1]);

    const Real maxRadius = maxRadiusFactor * radius;

    if (fieldSolver == "gmg") {
      run<FieldSolverGMG>(exact, maxRadius, numWorst, boxTagger);
    }
#ifdef CH_USE_PETSC
    else if (fieldSolver == "amg") {
      run<FieldSolverAMG>(exact, maxRadius, numWorst, boxTagger);
    }
#endif
    else {
      MayDay::Abort("DielectricVoid -- 'field_solver' must be 'gmg', or 'amg' in a build with USE_PETSC=TRUE");
    }
  }

  ChomboDischarge::finalize();
}
