/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_SigmaSolver.cpp
  @brief  Implementation of CD_SigmaSolver.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_SigmaSolver.H>
#include <CD_DataOps.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

constexpr int SigmaSolver::m_comp;
constexpr int SigmaSolver::m_nComp;

SigmaSolver::SigmaSolver()
{
  CH_TIME("SigmaSolver::SigmaSolver");

  // Default settings.
  this->setVerbosity(-1);
  this->setPhase(phase::gas);
  this->setRealm(Realm::Primal);
}

SigmaSolver::~SigmaSolver() { CH_TIME("SigmaSolver::~SigmaSolver"); }

const std::string
SigmaSolver::getRealm() const
{
  CH_TIME("SigmaSolver::getRealm");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::getRealm" << endl;
  }

  return m_realm;
}

void
SigmaSolver::setRealm(const std::string a_realm)
{
  CH_TIME("SigmaSolver::setRealm");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::setRealm" << endl;
  }

  m_realm = a_realm;
}

void
SigmaSolver::allocateInternals()
{
  CH_TIME("SigmaSolver::allocateInternals");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::allocateInternals" << endl;
  }

  m_amr->allocate(m_phi, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_flux, m_realm, m_phase, m_nComp);
}

void
SigmaSolver::preRegrid(const int a_lbase, const int a_oldFinestLevel)
{
  CH_TIME("SigmaSolver::preRegrid");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::preRegrid" << endl;
  }

  m_amr->allocate(m_cache, m_realm, m_phase, m_nComp);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    m_phi[lvl]->localCopyTo(*m_cache[lvl]);
  }
}

void
SigmaSolver::computeRHS(EBAMRIVData& a_rhs)
{
  CH_TIME("SigmaSolver::computeRHS");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::computeRHS" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    m_flux[lvl]->localCopyTo(*a_rhs[lvl]);
  }
}

void
SigmaSolver::deallocateInternals()
{
  CH_TIME("SigmaSolver::deallocateInternals");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::deallocateInternals" << endl;
  }

  m_amr->deallocate(m_phi);
  m_amr->deallocate(m_flux);
}

void
SigmaSolver::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("SigmaSolver::regrid");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::regrid" << endl;
  }

  // TLDR: This is pretty old code and I'm reluctant to change it even though it has a potential performance bottleneck here (due to the EBISLayout filling). What
  //       this code does is that we are essentially figuring out how to interpolate the data if a grid is added. Due to locality constraint, we have to coarsen
  //       the fine grids and copy the previous m_phi data onto those grids (the coarsened fine grids). This provides locality between the coarse grid and the fine grid
  //       and we can interpolate from that.

  const RefCountedPtr<EBIndexSpace> ebis = m_multifluidIndexSpace->getEBIndexSpace(m_phase);

  const int      ebghost = 4; // m_amr->getNumberOfEbGhostCells();
  const Interval interv(m_comp, m_comp);

  this->allocateInternals();

  DataOps::setValue(m_phi, 0.0);

  // These levels have never changed but the ownership might have, so we have to copy in order to ensure that data
  // ends up in the right place.
  for (int lvl = 0; lvl <= std::max(0, a_lmin - 1); lvl++) {
    m_cache[lvl]->copyTo(*m_phi[lvl]);
  }

  // These levels have changed and so we need to do the interpolation here.
  for (int lvl = std::max(1, a_lmin); lvl <= a_newFinestLevel; lvl++) {
    const DisjointBoxLayout& fineGrid   = m_amr->getGrids(m_realm)[lvl];
    const ProblemDomain&     coarDomain = m_amr->getDomains()[lvl - 1];
    const EBISLayout&        fineEBISL  = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const int                nref       = m_amr->getRefinementRatios()[lvl - 1];

    // Fill a coarsened grid and a coarsened ebisl
    EBISLayout        coarEBISL;
    DisjointBoxLayout coarGrid;
    coarsen(coarGrid, fineGrid, nref);
    ebis->fillEBISLayout(coarEBISL, coarGrid, coarDomain, ebghost);

    // Need extra storage for coarse data
    LayoutData<IntVectSet> sets(coarGrid);
    for (DataIterator dit = coarGrid.dataIterator(); dit.ok(); ++dit) {
      Box box = coarGrid.get(dit());
      box.grow(3);
      box &= coarDomain;

      const EBISBox&   ebisbox = coarEBISL[dit()];
      const IntVectSet ivs     = ebisbox.getIrregIVS(box);

      sets[dit()] = ivs;
    }

    // Allocate storage for coarsened data
    BaseIVFactory<Real>        ivfact(coarEBISL, sets);
    LevelData<BaseIVFAB<Real>> coFiData(coarGrid, m_nComp, m_phi[0]->ghostVect(), ivfact);

    // Copy the coarse data to the local view of the current fine data (we have already done the coarse level).
    EBLevelDataOps::setVal(coFiData, 0.0);
    m_phi[lvl - 1]->copyTo(coFiData);

    // Loop through coarse grid and interpolate to fine grid
    for (DataIterator dit = coarGrid.dataIterator(); dit.ok(); ++dit) {
      BaseIVFAB<Real>&       fineState   = (*m_phi[lvl])[dit()];
      const BaseIVFAB<Real>& coarState   = coFiData[dit()];
      const EBISBox&         coarEBISBox = coarEBISL[dit()];
      const EBISBox&         fineEBISBox = fineEBISL[dit()];

      const IntVectSet ivs     = fineState.getIVS();
      const EBGraph&   ebgraph = fineState.getEBGraph();

      // Iterate space for kernel.
      VoFIterator vofit(ivs, ebgraph);

      auto kernel = [&](const VolIndex& fineVof) -> void {
        const VolIndex coarVof  = fineEBISL.coarsen(fineVof, nref, dit());
        const Real     coarArea = coarEBISBox.bndryArea(coarVof);

        // Same sigma for all refined cells such that charge is conserved. This is equivalent to constant interpolation.
        const Vector<VolIndex> refinedVofs = coarEBISL.refine(coarVof, nref, dit());
        Real                   fineArea    = 0.0;
        for (int i = 0; i < refinedVofs.size(); i++) {
          fineArea += fineEBISBox.bndryArea(refinedVofs[i]);
        }

        // Initialize conserved charge. If there is no real area just set the state to zero.
        if (fineArea > 0.0 && coarArea > 0.0) {
          fineState(fineVof, m_comp) = coarState(coarVof, m_comp) * nref * coarArea / fineArea;
        }
        else {
          fineState(fineVof, m_comp) = 0.0;
        }
      };

      BoxLoops::loop(vofit, kernel);
    }

    // If data already exists, it takes precedence
    if (lvl <= a_oldFinestLevel) {
      m_cache[lvl]->copyTo(*m_phi[lvl]);
    }
  }

  this->resetCells(m_phi);
}

void
SigmaSolver::registerOperators()
{
  CH_TIME("SigmaSolver::registerOperators");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::registerOperators" << endl;
  }

  if (m_amr.isNull()) {
    MayDay::Error("SigmaSolver::registerOperators - need to set AmrMesh!");
  }
  else {
    m_amr->registerOperator(s_eb_coar_ave, m_realm, m_phase);
  }
}

void
SigmaSolver::resetCells(EBAMRIVData& a_data)
{
  CH_TIME("SigmaSolver::resetCells");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::resetCells" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[lvl];
    const MFLevelGrid&       mflg = *m_amr->getMFLevelGrid(m_realm)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      const Box        box  = dbl[dit()];
      BaseIVFAB<Real>& data = (*a_data[lvl])[dit()];

      IntVectSet     ivs     = data.getIVS();
      const EBGraph& ebgraph = data.getEBGraph();

      ivs -= mflg.interfaceRegion(box, dit());

      // Iteration space
      VoFIterator vofit(ivs, ebgraph);

      // Kernel
      auto kernel = [&](const VolIndex& vof) -> void {
        for (int comp = 0; comp < data.nComp(); comp++) {
          data(vofit(), comp) = 0.0;
        }
      };

      BoxLoops::loop(vofit, kernel);
    }
  }
}

void
SigmaSolver::setAmr(const RefCountedPtr<AmrMesh>& a_amr)
{
  CH_TIME("SigmaSolver::setAmr");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::setAmr" << endl;
  }

  m_amr = a_amr;
}

void
SigmaSolver::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry)
{
  CH_TIME("SigmaSolver::setComputationalGeometry");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::setComputationalGeometry" << endl;
  }

  m_computationalGeometry = a_computationalGeometry;
  m_multifluidIndexSpace  = m_computationalGeometry->getMfIndexSpace();
}

void
SigmaSolver::setPhase(phase::which_phase a_phase)
{
  CH_TIME("SigmaSolver::setPhase");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::setPhase" << endl;
  }

  m_phase = a_phase;
}

void
SigmaSolver::setSigma(const EBAMRIVData& a_sigma)
{
  CH_TIME("SigmaSolver::setSigma(ebamrivdata)");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::setSigma(ebamrivdata)" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    a_sigma[lvl]->localCopyTo(*m_phi[lvl]);
  }

  this->resetCells(m_phi);
}

void
SigmaSolver::setSigma(const Real a_sigma)
{
  CH_TIME("SigmaSolver::setSigma(Real)");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::setSigma(Real)" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    DataOps::setValue(*m_phi[lvl], a_sigma);
  }

  this->resetCells(m_phi);
}

void
SigmaSolver::setVerbosity(const int a_verbosity)
{
  CH_TIME("SigmaSolver::setVerbosity");

  m_verbosity = a_verbosity;

  if (m_verbosity > 5) {
    pout() << "SigmaSolver::setVerbosity" << endl;
  }
}

void
SigmaSolver::setTime(const int a_step, const Real a_time, const Real a_dt)
{
  CH_TIME("SigmaSolver::setTime");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::setTime" << endl;
  }

  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;
}

#ifdef CH_USE_HDF5
void
SigmaSolver::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const
{
  CH_TIME("SigmaSolver::writeCheckpointLevel");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::writeCheckpointLevel" << endl;
  }

  // Allocate some transient storage.
  EBCellFactory        fact(m_amr->getEBISLayout(m_realm, phase::gas)[a_level]);
  LevelData<EBCellFAB> scratch(m_amr->getGrids(m_realm)[a_level], 1, 3 * IntVect::Unit, fact);

  // Copy data into transient storage.
  DataOps::setValue(scratch, 0.0);
  DataOps::incr(scratch, *m_phi[a_level], 1.0);

  // Write state vector
  write(a_handle, scratch, "sigma");
}
#endif

#ifdef CH_USE_HDF5
void
SigmaSolver::readCheckpointLevel(HDF5Handle& a_handle, const int a_level)
{
  CH_TIME("SigmaSolver::readCheckpointLevel");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::readCheckpointLevel" << endl;
  }

  // TLDR: Chombo wants us to read into EBCellFABs so we allocate LevelData<EBCellFAB> and read into those. After that, we copy
  //       the result to our internal EBAMRIVFAB

  const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, phase::gas)[a_level];
  const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[a_level];

  // Transient that we read into.
  EBCellFactory        fact(ebisl);
  LevelData<EBCellFAB> scratch(dbl, 1, 3 * IntVect::Unit, fact);
  DataOps::setValue(scratch, 0.0);
  read<EBCellFAB>(a_handle, scratch, "sigma", dbl, Interval(0, 0), false);

  // Copy data to our EBAMRIVData
  DataOps::setValue(*m_phi[a_level], 0.0);
  DataOps::incr(*m_phi[a_level], scratch, 1.0);
}
#endif

void
SigmaSolver::writePlotData(EBAMRCellData& a_output, int& a_comp)
{
  CH_TIME("SigmaSolver::writePlotData");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::writePlotData" << endl;
  }

  // TLDR: We can't plot or write an EBAMRIVData data structure, so we copy our results to an EBAMRCellData so we can copy into a_output.

  // Allocate scratch storage that we can work with when we copy over to an EBAMRCellData.
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, 1);

  // Write sigma. Just copy over to scratch, and then copy over to a_output after that.
  DataOps::setValue(scratch, 0.0);
  DataOps::incr(scratch, m_phi, 1.0);
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    if (m_realm == a_output.getRealm()) {
      scratch[lvl]->localCopyTo(Interval(0, 0), *a_output[lvl], Interval(a_comp, a_comp));
    }
    else {
      scratch[lvl]->copyTo(Interval(0, 0), *a_output[lvl], Interval(a_comp, a_comp));
    }
  }
  a_comp++;

  // Write flux. Just the same thing as the code above.
  DataOps::setValue(scratch, 0.0);
  DataOps::incr(scratch, m_flux, 1.0);
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    if (m_realm == a_output.getRealm()) {
      scratch[lvl]->localCopyTo(Interval(0, 0), *a_output[lvl], Interval(a_comp, a_comp));
    }
    else {
      scratch[lvl]->copyTo(Interval(0, 0), *a_output[lvl], Interval(a_comp, a_comp));
    }
  }
  a_comp++;
}

int
SigmaSolver::getNumberOfPlotVariables()
{
  CH_TIME("SigmaSolver::getNumberOfPlotVariables");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::getNumberOfPlotVariables" << endl;
  }

  return 2;
}

Vector<std::string>
SigmaSolver::getPlotVariableNames() const
{
  CH_TIME("SigmaSolver::getNumberOfPlotVariables");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::getNumberOfPlotVariables" << endl;
  }

  Vector<std::string> ret(2);

  ret[0] = "surface charge density";
  ret[1] = "surface charge flux";

  return ret;
}

Real
SigmaSolver::computeCharge()
{
  CH_TIME("SigmaSolver::computeCharge");
  if (m_verbosity > 5) {
    pout() << "SigmaSolver::computeCharge" << endl;
  }

  // TLDR: Since we can coarsen conservatively, we compute coarsen the solution
  //       onto the coarsest level and do the integration there.

  m_amr->conservativeAverage(m_phi, m_realm, m_phase);

  Real charge = 0.0;

  const int                comp  = 0;
  const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[0];
  const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[0];
  const Real               dx    = m_amr->getDx()[0];

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box      box     = dbl[dit()];
    const EBISBox& ebisbox = ebisl[dit()];

    const EBGraph&   ebgraph  = ebisbox.getEBGraph();
    const IntVectSet irregIVS = ebisbox.getIrregIVS(box);

    // Kernel space
    VoFIterator vofit(irregIVS, ebgraph);

    // Kernel
    auto kernel = [&](const VolIndex& vof) -> void {
      const Real& area = ebisbox.bndryArea(vof);

      charge += area * (*m_phi[0])[dit()](vof, comp);
    };

    BoxLoops::loop(vofit, kernel);
  }

  DataOps::sum(charge); // Parallell sum.

  // Return charge
  return charge * std::pow(dx, SpaceDim - 1);
}

EBAMRIVData&
SigmaSolver::getPhi()
{
  return m_phi;
}

EBAMRIVData&
SigmaSolver::getFlux()
{
  return m_flux;
}

#include <CD_NamespaceFooter.H>
