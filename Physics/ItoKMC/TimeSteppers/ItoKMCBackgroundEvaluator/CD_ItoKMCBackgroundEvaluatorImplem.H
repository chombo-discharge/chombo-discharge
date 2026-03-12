/* chombo-discharge
 * Copyright © 2025 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCBackgroundEvaluatorImplem.H
  @brief  Implementation of CD_ItoKMCBackgroundEvaluator.H
  @author Robert Marskar
*/

#ifndef CD_ItoKMCBackgroundEvaluatorImplem_H
#define CD_ItoKMCBackgroundEvaluatorImplem_H

// Our includes
#include <CD_ItoKMCBackgroundEvaluator.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

template <typename I, typename C, typename R, typename F>
ItoKMCBackgroundEvaluator<I, C, R, F>::ItoKMCBackgroundEvaluator(RefCountedPtr<ItoKMCPhysics>& a_physics) noexcept
  : ItoKMCGodunovStepper<I, C, R, F>(a_physics, false)
{
  CH_TIME("ItoKMCBackgroundEvaluator::ItoKMCBackgroundEvaluator");

  // Set default values
  m_maxFieldExitCrit   = std::numeric_limits<Real>::max();
  m_relFieldExitCrit   = std::numeric_limits<Real>::max();
  m_maxFieldChange     = -1.0;
  m_relFieldChange     = -1.0;
  m_electrodeCharge    = 0.0;
  m_ohmicCharge        = 0.0;
  m_opticalSolver      = "2PN2";
  m_maxInitialTimeStep = 0.0;

  // Override the name and re-parse ALL options (ItoKMCStepper + GodunovStepper) using the new prefix
  this->m_name = "ItoKMCBackgroundEvaluator";
  this->parseOptions();

  // Parse BackgroundEvaluator-specific options
  ParmParse pp(this->m_name.c_str());
  pp.query("max_field_exit_crit", m_maxFieldExitCrit);
  pp.query("rel_field_exit_crit", m_relFieldExitCrit);
  pp.query("optical_solver", m_opticalSolver);
  pp.query("max_initial_time_step", m_maxInitialTimeStep);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCBackgroundEvaluator<I, C, R, F>::postInitialize() noexcept
{
  CH_TIME("ItoKMCBackgroundEvaluator::postInitialize");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCBackgroundEvaluator::postInitialize" << endl;
  }

  ItoKMCGodunovStepper<I, C, R, F>::postInitialize();

  this->computeBackgroundField();
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCBackgroundEvaluator<I, C, R, F>::postRegrid() noexcept
{
  CH_TIME("ItoKMCBackgroundEvaluator::postRegrid");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCBackgroundEvaluator::postRegrid" << endl;
  }

  ItoKMCGodunovStepper<I, C, R, F>::postRegrid();

  this->computeBackgroundField();
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCBackgroundEvaluator<I, C, R, F>::postCheckpointSetup() noexcept
{
  CH_TIME("ItoKMCBackgroundEvaluator::postCheckpointSetup");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCBackgroundEvaluator::postCheckpointSetup" << endl;
  }

  ItoKMCGodunovStepper<I, C, R, F>::postCheckpointSetup();

  this->computeBackgroundField();
}

template <typename I, typename C, typename R, typename F>
Real
ItoKMCBackgroundEvaluator<I, C, R, F>::advance(const Real a_dt)
{
  CH_TIME("ItoKMCBackgroundEvaluator::advance");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCBackgroundEvaluator::advance" << endl;
  }

  const Real actualDt = ItoKMCGodunovStepper<I, C, R, F>::advance(a_dt);

  const std::pair<Real, Real> fieldChange      = this->evaluateSpaceChargeEffects();
  const Real                  electrodeQ       = this->integrateElectrodeSurfaceCharge();
  const Real                  ohmicQ           = this->integrateOhmicCharge();
  const std::pair<Real, Real> opticalEmissions = this->integrateOpticalExcitations();

  m_relFieldChange  = fieldChange.first;
  m_maxFieldChange  = fieldChange.second;
  m_electrodeCharge = electrodeQ;
  m_ohmicCharge     = ohmicQ;
  m_sumOpticalPhi   = opticalEmissions.first;
  m_sumOpticalSrc   = opticalEmissions.second;

  bool abortCondition = false;
  if ((m_maxFieldChange > m_maxFieldExitCrit) && (m_maxFieldExitCrit > 0.0)) {
    pout() << "ItoKMCBackgroundEvaluator -- stopping because the max field (in any cell) "
           << "changed by more than specified threshold (max_field_exit_crit)" << endl;
    abortCondition = true;
  }
  if ((m_relFieldChange > m_relFieldExitCrit) && (m_relFieldExitCrit > 0.0)) {
    pout() << "ItoKMCBackgroundEvaluator -- stopping because the relative field change in a cell "
           << "is larger than specified threshold (rel_field_exit_crit)" << endl;
    abortCondition = true;
  }
  if (abortCondition) {
    this->m_keepGoing = false;
  }

  return actualDt;
}

template <typename I, typename C, typename R, typename F>
Real
ItoKMCBackgroundEvaluator<I, C, R, F>::computeDt()
{
  CH_TIME("ItoKMCBackgroundEvaluator::computeDt");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCBackgroundEvaluator::computeDt" << endl;
  }

  Real newDt = ItoKMCGodunovStepper<I, C, R, F>::computeDt();

  if (this->m_timeStep == 0 && m_maxInitialTimeStep > 0.0 && newDt > m_maxInitialTimeStep)
    newDt = m_maxInitialTimeStep;

  return newDt;
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCBackgroundEvaluator<I, C, R, F>::printStepReport() noexcept
{
  CH_TIME("ItoKMCBackgroundEvaluator::printStepReport");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCBackgroundEvaluator::printStepReport" << endl;
  }

  std::ios::fmtflags oldFlags = pout().flags();

  ItoKMCGodunovStepper<I, C, R, F>::printStepReport();

  // Append the step report.
  const std::string whitespace = "                                   ";
  pout().precision(12);
  pout() << whitespace + "Delta E(max)      = " << 100.0 * m_maxFieldChange << " (%)" << endl;
  pout() << whitespace + "Delta E(rel)      = " << 100.0 * m_relFieldChange << " (%)" << endl;
  pout() << whitespace + "Q (electrode)     = " << std::scientific << m_electrodeCharge << " (C)" << endl;
  pout() << whitespace + "Q (ohmic)         = " << std::scientific << m_ohmicCharge << " (C)" << endl;
  pout() << whitespace + "Sum (phi_optical) = " << std::scientific << m_sumOpticalPhi << endl;
  pout() << whitespace + "Sum (src_optical) = " << std::scientific << m_sumOpticalSrc << " (1/s)" << endl;

  pout().flags(oldFlags);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCBackgroundEvaluator<I, C, R, F>::computeBackgroundField() noexcept
{
  CH_TIME("ItoKMCBackgroundEvaluator::computeBackgroundField");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCBackgroundEvaluator::computeBackgroundField" << endl;
  }

  // Allocate the background field, storage for the potential, storage for the space-charge, and storage for the
  // surface charge.
  MFAMRCellData phi;
  MFAMRCellData rho;
  EBAMRIVData   sigma;

  (this->m_amr)->allocate(phi, this->m_fluidRealm, 1);
  (this->m_amr)->allocate(rho, this->m_fluidRealm, 1);
  (this->m_amr)->allocate(sigma, this->m_fluidRealm, this->m_plasmaPhase, 1);
  (this->m_amr)->allocate(m_backgroundField, this->m_fluidRealm, SpaceDim);

  // Copy the potential over from the field solver so we have a better initial guess for the potential. The
  // space-charge and surface charge must be zero.
  (this->m_amr)->copyData(phi, (this->m_fieldSolver)->getPotential());

  DataOps::setValue(rho, 0.0);
  DataOps::setValue(sigma, 0.0);

  // Reset the field solver permittivities. Note that this discards the semi-implicit coefficients for the field update,
  // but these are refilled on the next solver time step anyways.
  (this->m_fieldSolver)->setPermittivities();

  const MFAMRCellData& permCell = (this->m_fieldSolver)->getPermittivityCell();
  const MFAMRFluxData& permFace = (this->m_fieldSolver)->getPermittivityFace();
  const MFAMRIVData&   permEB   = (this->m_fieldSolver)->getPermittivityEB();

  (this->m_fieldSolver)->setSolverPermittivities(permCell, permFace, permEB);

  // Solve the damn thing, and then compute the electric field onto m_backgroundField.
  (this->m_fieldSolver)->solve(phi, rho, sigma, false);
  (this->m_fieldSolver)->computeElectricField(m_backgroundField, phi);
}

template <typename I, typename C, typename R, typename F>
std::pair<Real, Real>
ItoKMCBackgroundEvaluator<I, C, R, F>::evaluateSpaceChargeEffects() noexcept
{
  CH_TIME("ItoKMCBackgroundEvaluator::evaluateSpaceChargeEffects");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCBackgroundEvaluator::evaluateSpaceChargeEffects" << endl;
  }

  EBAMRCellData poissonNorm;
  EBAMRCellData laplaceNorm;
  EBAMRCellData poissonField;
  EBAMRCellData laplaceField;
  EBAMRCellData deltaE;

  (this->m_amr)->allocate(poissonNorm, this->m_fluidRealm, this->m_plasmaPhase, 1);
  (this->m_amr)->allocate(laplaceNorm, this->m_fluidRealm, this->m_plasmaPhase, 1);
  (this->m_amr)->allocate(poissonField, this->m_fluidRealm, this->m_plasmaPhase, SpaceDim);
  (this->m_amr)->allocate(laplaceField, this->m_fluidRealm, this->m_plasmaPhase, SpaceDim);
  (this->m_amr)->allocate(deltaE, this->m_fluidRealm, this->m_plasmaPhase, 1);

  // Get the fields in the plasma phase -- the actual field and the background field.
  const MFAMRCellData electricField = (this->m_fieldSolver)->getElectricField();

  DataOps::copy(poissonField, (this->m_amr)->alias(this->m_plasmaPhase, electricField));
  DataOps::copy(laplaceField, (this->m_amr)->alias(this->m_plasmaPhase, m_backgroundField));

  (this->m_amr)->interpToCentroids(poissonField, this->m_fluidRealm, this->m_plasmaPhase);
  (this->m_amr)->interpToCentroids(laplaceField, this->m_fluidRealm, this->m_plasmaPhase);

  // Compute the field magnitude, and then the relative change in field magnitude.
  DataOps::vectorLength(poissonNorm, poissonField);
  DataOps::vectorLength(laplaceNorm, laplaceField);

  DataOps::copy(deltaE, poissonNorm);
  DataOps::incr(deltaE, laplaceNorm, -1.0);
  DataOps::divideFallback(deltaE, laplaceNorm, 0.0);

  (this->m_amr)->arithmeticAverage(deltaE, this->m_fluidRealm, this->m_plasmaPhase);
  (this->m_amr)->interpGhost(deltaE, this->m_fluidRealm, this->m_plasmaPhase);

  // Iterate through the grid cells and figure out where the field change the most.
  //
  // PS: I'm doing this with a direct loop because DataOps::getMaxMin will do ALL cells, including ones that are
  //     covered by a finer level. But we only want to evaluate the valid region.
  Real relChange = -1.0;
  Real maxChange = -1.0;

  Real maxSpaceChargeField = 0.0;
  Real minSpaceChargeField = 0.0;
  Real maxBackgroundField  = 0.0;
  Real minBackgroundField  = 0.0;

  DataOps::getMaxMin(maxSpaceChargeField, minSpaceChargeField, poissonNorm, 0);
  DataOps::getMaxMin(maxBackgroundField, minBackgroundField, laplaceNorm, 0);

  for (int lvl = 0; lvl <= (this->m_amr)->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = (this->m_amr)->getGrids(this->m_fluidRealm)[lvl];
    const DataIterator&      dit   = dbl.dataIterator();
    const EBISLayout&        ebisl = (this->m_amr)->getEBISLayout(this->m_fluidRealm, this->m_plasmaPhase)[lvl];

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(max : relChange)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex&     din        = dit[mybox];
      const Box            box        = dbl[din];
      const EBISBox&       ebisbox    = ebisl[din];
      const BaseFab<bool>& validCells = (*(this->m_amr)->getValidCells(this->m_fluidRealm)[lvl])[din];

      const EBCellFAB& data    = (*deltaE[lvl])[din];
      const FArrayBox& dataReg = data.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        if (validCells(iv) && ebisbox.isRegular(iv)) {
          relChange = std::max(relChange, std::abs(dataReg(iv, 0)));
        }
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        if (validCells(vof.gridIndex()) && ebisbox.isIrregular(vof.gridIndex())) {
          relChange = std::max(relChange, std::abs(data(vof, 0)));
        }
      };

      VoFIterator& vofit = (*(this->m_amr)->getVofIterator(this->m_fluidRealm, this->m_plasmaPhase)[lvl])[din];

      BoxLoops::loop(box, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }

  relChange = ParallelOps::max(relChange);
  maxChange = maxSpaceChargeField / maxBackgroundField - 1.0;

  return std::make_pair(relChange, maxChange);
}

template <typename I, typename C, typename R, typename F>
Real
ItoKMCBackgroundEvaluator<I, C, R, F>::integrateElectrodeSurfaceCharge() const noexcept
{
  CH_TIME("ItoKMCBackgroundEvaluator::integrateElectrodeSurfaceCharge");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCBackgroundEvaluator::integrateElectrodeSurfaceCharge" << endl;
  }

  Real Enorm = 0.0;

  const Real eps = Units::eps0 * this->m_computationalGeometry->getGasPermittivity();

  const MFAMRCellData& electricFieldCell = (this->m_fieldSolver)->getElectricField();
  const MFAMRCellData& permittivityCell  = (this->m_fieldSolver)->getPermittivityCell();
  const EBAMRCellData  permittivity      = (this->m_amr)->alias(this->m_plasmaPhase, permittivityCell);

  // Allocates data and computes the electric field on the cell centroid(s).
  EBAMRCellData electricField;
  (this->m_amr)->allocate(electricField, this->m_fluidRealm, this->m_plasmaPhase, SpaceDim);
  DataOps::copy(electricField, (this->m_amr)->alias(this->m_plasmaPhase, electricFieldCell));
  (this->m_amr)->interpToCentroids(electricField, this->m_fluidRealm, this->m_plasmaPhase);

  for (int lvl = 0; lvl <= (this->m_amr)->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = (this->m_amr)->getGrids(this->m_fluidRealm)[lvl];
    const EBISLayout&        ebisl = (this->m_amr)->getEBISLayout(this->m_fluidRealm, this->m_plasmaPhase)[lvl];
    const DataIterator&      dit   = dbl.dataIterator();
    const Real               dx    = (this->m_amr)->getDx()[lvl];
    const Real               dA    = std::pow(dx, SpaceDim - 1);

    const int numBoxes = dit.size();

#pragma omp parallel for schedule(runtime) reduction(+ : Enorm)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex&     din        = dit[mybox];
      const EBISBox&       ebisbox    = ebisl[din];
      const EBCellFAB&     E          = (*electricField[lvl])[din];
      const EBCellFAB&     eps        = (*permittivity[lvl])[din];
      const BaseFab<bool>& validCells = (*(this->m_amr)->getValidCells(this->m_fluidRealm)[lvl])[din];

      auto kernel = [&](const VolIndex& vof) -> void {
        if (ebisbox.isIrregular(vof.gridIndex()) && validCells(vof.gridIndex(), 0)) {
          const RealVect normalVector = ebisbox.normal(vof);
          const RealVect e            = RealVect(D_DECL(E(vof, 0), E(vof, 1), E(vof, 2)));
          const Real     areaFrac     = ebisbox.bndryArea(vof);

          Enorm += eps(vof, 0) * e.dotProduct(normalVector) * areaFrac * dA;
        }
      };

      VoFIterator& vofit = (*(this->m_amr)->getVofIterator(this->m_fluidRealm, this->m_plasmaPhase)[lvl])[din];

      BoxLoops::loop(vofit, kernel);
    }
  }

  return eps * ParallelOps::sum(Enorm);
}

template <typename I, typename C, typename R, typename F>
Real
ItoKMCBackgroundEvaluator<I, C, R, F>::integrateOhmicCharge() const noexcept
{
  CH_TIME("ItoKMCBackgroundEvaluator::integrateOhmicCharge");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCBackgroundEvaluator::integrateOhmicCharge" << endl;
  }

  return 0.0;
}

template <typename I, typename C, typename R, typename F>
std::pair<Real, Real>
ItoKMCBackgroundEvaluator<I, C, R, F>::integrateOpticalExcitations() const noexcept
{
  CH_TIME("ItoKMCBackgroundEvaluator::integrateOpticalExcitations");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCBackgroundEvaluator::integrateOpticalExcitations" << endl;
  }

  Real sumPhi = 0.0;
  Real sumSrc = 0.0;

  for (CdrIterator<CdrSolver> solverIt = (this->m_cdr)->iterator(); solverIt.ok(); ++solverIt) {
    const RefCountedPtr<CdrSolver>& solver = solverIt();

    if (solver->getName() == m_opticalSolver) {
      for (int lvl = 0; lvl <= (this->m_amr)->getFinestLevel(); lvl++) {
        const DisjointBoxLayout& dbl   = (this->m_amr)->getGrids(this->m_fluidRealm)[lvl];
        const EBISLayout&        ebisl = (this->m_amr)->getEBISLayout(this->m_fluidRealm, this->m_plasmaPhase)[lvl];
        const DataIterator&      dit   = dbl.dataIterator();
        const Real               dx    = (this->m_amr)->getDx()[lvl];
        const Real               dV    = std::pow(dx, SpaceDim);

        const int numBoxes = dit.size();

#pragma omp parallel for schedule(runtime) reduction(+ : sumPhi, sumSrc)
        for (int mybox = 0; mybox < numBoxes; mybox++) {
          const DataIndex&     din        = dit[mybox];
          const Box            cellBox    = dbl[din];
          const EBISBox&       ebisbox    = ebisl[din];
          const BaseFab<bool>& validCells = (*(this->m_amr)->getValidCells(this->m_fluidRealm)[lvl])[din];

          const EBCellFAB& phi = (*(solver->getPhi()[lvl]))[din];
          const EBCellFAB& src = (*(solver->getSource()[lvl]))[din];

          const FArrayBox& phiReg = phi.getFArrayBox();
          const FArrayBox& srcReg = src.getFArrayBox();

          auto regularKernel = [&](const IntVect& iv) -> void {
            if (ebisbox.isRegular(iv) && validCells(iv)) {
              sumPhi += phiReg(iv, 0) * dV;
              sumSrc += srcReg(iv, 0) * dV;
            }
          };

          auto irregularKernel = [&](const VolIndex& vof) -> void {
            if (ebisbox.isIrregular(vof.gridIndex()) && validCells(vof.gridIndex(), 0)) {
              const Real kappa = ebisbox.volFrac(vof);

              sumPhi += phi(vof, 0) * dV;
              sumSrc += src(vof, 0) * dV;
            }
          };

          VoFIterator& vofit = (*(this->m_amr)->getVofIterator(this->m_fluidRealm, this->m_plasmaPhase)[lvl])[din];

          BoxLoops::loop(cellBox, regularKernel);
          BoxLoops::loop(vofit, irregularKernel);
        }
      }
    }
  }

  sumPhi = ParallelOps::sum(sumPhi);
  sumSrc = ParallelOps::sum(sumSrc);

  return std::make_pair(sumPhi, sumSrc);
}

#include <CD_NamespaceFooter.H>

#endif
