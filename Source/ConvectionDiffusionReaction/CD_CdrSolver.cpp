/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrSolver.cpp
  @brief  Implementation of CD_CdrSolver.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <EBAMRIO.H>
#include <EBArith.H>

// Our includes
#include <CD_CdrSolver.H>
#include <CD_DataOps.H>
#include <CD_BoxLoops.H>
#include <CD_ParallelOps.H>
#include <CD_DischargeIO.H>
#include <CD_Random.H>
#include <CD_NamespaceHeader.H>

constexpr int CdrSolver::m_comp;
constexpr int CdrSolver::m_nComp;

CdrSolver::CdrSolver()
{

  // Default options.
  m_verbosity    = -1;
  m_name         = "CdrSolver";
  m_className    = "CdrSolver";
  m_regridSlopes = true;

  this->setRealm(Realm::Primal);
  this->setDefaultDomainBC(); // Set default domain BCs (wall)
}

CdrSolver::~CdrSolver()
{}

RefCountedPtr<CdrSpecies>&
CdrSolver::getSpecies() noexcept
{
  CH_TIME("CdrSolver::getSpecies()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getSpecies()" << endl;
  }

  return m_species;
}

const RefCountedPtr<CdrSpecies>&
CdrSolver::getSpecies() const noexcept
{
  CH_TIME("CdrSolver::getSpecies()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getSpecies()" << endl;
  }

  return m_species;
}

void
CdrSolver::setDefaultDomainBC()
{
  CH_TIME("CdrSolver::setDefaultDomainBC()");
  if (m_verbosity > 5) {
    pout() << m_name + "::setDefaultDomainBC()" << endl;
  }

  // TLDR: This sets the domain boundary condition to be a wall BC (no incoming/outgoing mass).

  // Lambda function for wall bc -- mostly left in place so I can remind myself how to do this.
  auto zero = [](const RealVect a_position, const Real a_time) -> Real {
    return 0.0;
  };

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      const CdrDomainBC::DomainSide domainSide = std::make_pair(dir, sit());

      this->setDomainBcType(domainSide, CdrDomainBC::BcType::Wall);
      this->setDomainBcFunction(domainSide, zero);
    }
  }
}

void
CdrSolver::setDomainBcType(const CdrDomainBC::DomainSide a_domainSide, const CdrDomainBC::BcType a_bcType)
{
  CH_TIME("CdrSolver::setDomainBcType(CdrDomainBC::DomainSide, CdrDomainBC::BcType)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setDomainBcType(CdrDomainBC::DomainSide, CdrDomainBC::BcType)" << endl;
  }

  m_domainBC.setBcType(a_domainSide, a_bcType);
}

void
CdrSolver::setDomainBcFunction(const CdrDomainBC::DomainSide   a_domainSide,
                               const CdrDomainBC::FluxFunction a_fluxFunction)
{
  CH_TIME("CdrSolver::setDomainBcFunction(CdrDomainBC::DomainSide, CdrDomainBC::FluxFunction)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setDomainBcFunction(CdrDomainBC::DomainSide, CdrDomainBC::FluxFunction)" << endl;
  }

  m_domainBC.setBcFunction(a_domainSide, a_fluxFunction);
}

std::string
CdrSolver::getName() const
{
  CH_TIME("CdrSolver::getName()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getName()" << endl;
  }

  return m_name;
}

std::string
CdrSolver::getRealm() const
{
  CH_TIME("CdrSolver::getRealm()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getRealm()" << endl;
  }

  return m_realm;
}

void
CdrSolver::setRealm(const std::string a_realm)
{
  CH_TIME("CdrSolver::setRealm(std::string)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setRealm(std::string)" << endl;
  }

  m_realm.assign(a_realm);
}

Vector<std::string>
CdrSolver::getPlotVariableNames() const
{
  CH_TIME("CdrSolver::getPlotVariableNames()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getPlotVariableNames()" << endl;
  }

  // TLDR: Possible plot variables is the density (m_phi), diffusion coefficient, source term, velocity, and eb flux. This
  //       function returns the associated plot variable names, and will be used in the plot files.

  Vector<std::string> plotVarNames(0);

  if (m_plotPhi) {
    plotVarNames.push_back(m_name + " phi");
  }
  if (m_plotDiffusionCoefficient && m_isDiffusive) {
    plotVarNames.push_back(m_name + " diffusion_coefficient");
  }
  if (m_plotSource) {
    plotVarNames.push_back(m_name + " source");
  }
  if (m_plotVelocity && m_isMobile) {
    plotVarNames.push_back("x-Velocity " + m_name);
  }
  if (m_plotVelocity && m_isMobile) {
    plotVarNames.push_back("y-Velocity " + m_name);
  }
  if (m_plotVelocity && m_isMobile && SpaceDim == 3) {
    plotVarNames.push_back("z-Velocity " + m_name);
  }
  if (m_plotEbFlux && (m_isMobile || m_isDiffusive)) {
    plotVarNames.push_back(m_name + " eb_flux");
  }

  return plotVarNames;
}

int
CdrSolver::getNumberOfPlotVariables() const
{
  CH_TIME("CdrSolver::getNumberOfPlotVariables()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getNumberOfPlotVariables()" << endl;
  }

  // TLDR: Possible plot variables is the density (m_phi), diffusion coefficient, source term, velocity, and eb flux. This
  //       function returns the number of variables that this sum up to (a vector field is 2/3 variables in 2D/3D)

  int numPlotVars = 0;

  if (m_plotPhi) {
    numPlotVars = numPlotVars + 1;
  }
  if (m_plotDiffusionCoefficient && m_isDiffusive) {
    numPlotVars = numPlotVars + 1;
  }
  if (m_plotSource) {
    numPlotVars = numPlotVars + 1;
  }
  if (m_plotVelocity && m_isMobile) {
    numPlotVars = numPlotVars + SpaceDim;
  }
  if (m_plotEbFlux && m_isMobile) {
    numPlotVars = numPlotVars + 1;
  }

  return numPlotVars;
}

void
CdrSolver::advanceEuler(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const Real a_dt)
{
  CH_TIME("CdrSolver::advanceEuler(EBAMRCellData, EBAMRCellData, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::advanceEuler(EBAMRCellData, EBAMRCellData, Real)" << endl;
  }

  CH_assert(a_newPhi[0]->nComp() == 1);
  CH_assert(a_oldPhi[0]->nComp() == 1);

  // TLDR: We are solving phi^(k+1) - phi^k - dt*Div(D*Grad(phi^(k+1))) = 0.0. We create a source term = 0 and call the implementation
  //       version.
  if (m_isDiffusive) {
    EBAMRCellData src;
    m_amr->allocate(src, m_realm, m_phase, m_nComp);
    DataOps::setValue(src, 0.0);

    this->advanceEuler(a_newPhi, a_oldPhi, src, a_dt);
  }
  else {
    DataOps::copy(a_newPhi, a_oldPhi);
  }
}

void
CdrSolver::advanceCrankNicholson(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const Real a_dt)
{
  CH_TIME("CdrSolver::advanceCrankNicholson(EBAMRCellData, EBAMRCellData, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::advanceCrankNicholson(EBAMRCellData, EBAMRCellData, Real)" << endl;
  }

  CH_assert(a_newPhi[0]->nComp() == 1);
  CH_assert(a_oldPhi[0]->nComp() == 1);

  // TLDR: We are solving phi^(k+1) - phi^k - dt*Div(D*Grad(phi^(k+1))) = 0.0. We create a source term = 0 and call the implementation
  //       version.
  if (m_isDiffusive) {
    EBAMRCellData src;
    m_amr->allocate(src, m_realm, m_phase, m_nComp);
    DataOps::setValue(src, 0.0);

    this->advanceCrankNicholson(a_newPhi, a_oldPhi, src, a_dt);
  }
  else {
    DataOps::copy(a_newPhi, a_oldPhi);
  }
}

void
CdrSolver::allocate()
{
  CH_TIME("CdrSolver::allocate()");
  if (m_verbosity > 5) {
    pout() << m_name + "::allocate()" << endl;
  }

  // TLDR: This allocates a number of fields for this solver. We wish to trim memory if we can, so we don't allocate
  //       storage for velocities/diffusion coefficients if those fields are not going to be used. Somewhat confusingly,
  //       we allocate pointers for these fields but no memory blocks (for interface reasons).

  // These three are allocated no matter what -- we will always need the state (and probably also the source term)
  m_amr->allocate(m_phi, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_source, m_realm, m_phase, m_nComp);

  DataOps::setValue(m_phi, 0.0);
  DataOps::setValue(m_source, 0.0);

  // Only allocate memory for cell-centered and face-centered velocities if the solver is mobile. Otherwise, allocate
  // a NULL pointer that we can pass around in TimeStepper in order to handle special cases
  if (m_isMobile) {
    m_amr->allocate(m_cellVelocity, m_realm, m_phase, SpaceDim);
    m_amr->allocate(m_faceVelocity, m_realm, m_phase, m_nComp);
    m_amr->allocate(m_faceStates, m_realm, m_phase, m_nComp);

    DataOps::setValue(m_cellVelocity, 0.0);
    DataOps::setValue(m_faceVelocity, 0.0);
    DataOps::setValue(m_faceStates, 0.0);
  }
  else {
    m_amr->allocatePointer(m_faceVelocity, m_realm);
    m_amr->allocatePointer(m_cellVelocity, m_realm);
    m_amr->allocatePointer(m_faceStates, m_realm);
  }

  // Only allocate memory for diffusion coefficients if we need it. Otherwise, allocate a NULL pointer that we can
  // pass around pointers in case we need it. This might seem confusing but its easier to resize those vectors
  // and set nullpointers when we iterate through grid levels and patches.
  if (m_isDiffusive) {
    m_amr->allocate(m_cellCenteredDiffusionCoefficient, m_realm, m_phase, m_nComp);
    m_amr->allocate(m_faceCenteredDiffusionCoefficient, m_realm, m_phase, m_nComp);
    m_amr->allocate(m_ebCenteredDiffusionCoefficient, m_realm, m_phase, m_nComp);

    DataOps::setValue(m_cellCenteredDiffusionCoefficient, 0.0);
    DataOps::setValue(m_faceCenteredDiffusionCoefficient, 0.0);
    DataOps::setValue(m_ebCenteredDiffusionCoefficient, 0.0);
  }
  else {
    m_amr->allocatePointer(m_cellCenteredDiffusionCoefficient, m_realm);
    m_amr->allocatePointer(m_faceCenteredDiffusionCoefficient, m_realm);
    m_amr->allocatePointer(m_ebCenteredDiffusionCoefficient, m_realm);
  }

  // Allocate stuff for holding fluxes -- this data is used when computing advection and diffusion fluxes.
  if (m_isDiffusive || m_isMobile) {
  }

  // These don't consume (much) memory so we always allocate them.
  m_amr->allocate(m_ebFlux, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_ebZero, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_domainFlux, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_massDifference, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_nonConservativeDivG, m_realm, m_phase, m_nComp);

  DataOps::setValue(m_ebFlux, 0.0);
  DataOps::setValue(m_ebZero, 0.0);
  DataOps::setValue(m_domainFlux, 0.0);
  DataOps::setValue(m_massDifference, 0.0);
  DataOps::setValue(m_nonConservativeDivG, 0.0);

  // Define interpolation stencils
  this->defineInterpolationStencils();
}

void
CdrSolver::deallocate()
{
  CH_TIME("CdrSolver::deallocate()");
  if (m_verbosity > 5) {
    pout() << m_name + "::deallocate()" << endl;
  }

  m_phi.clear();
  m_source.clear();
  m_faceVelocity.clear();
  m_faceStates.clear();
  m_cellVelocity.clear();
  m_ebFlux.clear();
  m_cellCenteredDiffusionCoefficient.clear();
  m_faceCenteredDiffusionCoefficient.clear();
  m_ebCenteredDiffusionCoefficient.clear();
}

void
CdrSolver::averageVelocityToFaces()
{
  CH_TIME("CdrSolver::averageVelocityToFaces()");
  if (m_verbosity > 5) {
    pout() << m_name + "::averageVelocityToFaces()" << endl;
  }

  this->averageVelocityToFaces(m_faceVelocity, m_cellVelocity); // Average velocities to face centers for all levels
}

void
CdrSolver::averageVelocityToFaces(EBAMRFluxData& a_faceVelocity, const EBAMRCellData& a_cellVelocity)
{
  CH_TIME("CdrSolver::averageVelocityToFaces(EBAMRFluxData, EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::averageVelocityToFaces(EBAMRFluxData, EBAMRCellData)" << endl;
  }

  CH_assert(a_faceVelocity[0]->nComp() == 1);
  CH_assert(a_cellVelocity[0]->nComp() == SpaceDim);

#ifndef NDEBUG
  // Put in something huge in debug code so we can catch if the code breaks.
  DataOps::setValue(a_faceVelocity, std::numeric_limits<Real>::max());
#endif

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const int tanGhosts = 1;

    DataOps::averageCellVelocityToFaceVelocity(*a_faceVelocity[lvl],
                                               *a_cellVelocity[lvl],
                                               m_amr->getDomains()[lvl],
                                               tanGhosts);
  }
}

void
CdrSolver::preRegrid(const int a_lmin, const int a_oldFinestLevel)
{
  CH_TIME("CdrSolver::preRegrid(int, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::preRegrid(int, int)" << endl;
  }

  m_amr->allocate(m_cachePhi, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_cacheSource, m_realm, m_phase, m_nComp);

  m_amr->copyData(m_cachePhi, m_phi);
  m_amr->copyData(m_cacheSource, m_source);

  this->deallocate();
}

void
CdrSolver::computeDivG(EBAMRCellData&     a_divG,
                       EBAMRFluxData&     a_G,
                       const EBAMRIVData& a_ebFlux,
                       const bool         a_conservativeOnly)
{
  CH_TIME("CdrSolver::computeDivG(EBAMRCellData, EBAMRFluxData, EBAMRIVData, bool)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDivG(EBAMRCellData, EBAMRFluxData, EBAMRIVData, bool)" << endl;
  }

  CH_assert(a_divG[0]->nComp() == 1);
  CH_assert(a_G[0]->nComp() == 1);
  CH_assert(a_ebFlux[0]->nComp() == 1);

  // TLDR: This routine computes a finite volume approximation to Div(G) where G is a flux, stored on face centers and eb faces. The routine uses
  //       flux matching on refinement boundaries and the so-called hybrid divergence in the cut-cells. The mass which is missed by the hybrid
  //       divergence is smooshed back in through through redistribution.

  DataOps::setValue(a_divG, 0.0);

  // Interpolate fluxes to centroids and then coarsen them so we have consistent fluxes.
  this->interpolateFluxToFaceCentroids(a_G);
  m_amr->conservativeAverage(a_G, m_realm, m_phase);

  // Make the conservative divergence without kappa-division.
  this->conservativeDivergenceNoKappaDivision(a_divG, a_G, a_ebFlux);

  // Compute hybrid divergence.
  if (!a_conservativeOnly) {
    // Compute the non-conservative divergence
    this->nonConservativeDivergence(m_nonConservativeDivG, a_divG);

    // a_divG becomes hybrid divergence. Mass diff computed.
    this->hybridDivergence(a_divG, m_massDifference, m_nonConservativeDivG);

    if (m_whichRedistribution != Redistribution::None) {

      Vector<RefCountedPtr<EBFluxRedistribution>>& redistOps = m_amr->getRedistributionOp(m_realm, m_phase);

      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const Real     scale     = 1.0;
        const Interval variables = Interval(0, 0);
        const bool     hasCoar   = lvl > 0;
        const bool     hasFine   = lvl < m_amr->getFinestLevel();

        if (hasCoar) {
          redistOps[lvl]->redistributeCoar(*a_divG[lvl - 1], *m_massDifference[lvl], scale, variables);
        }

        redistOps[lvl]->redistributeLevel(*a_divG[lvl], *m_massDifference[lvl], scale, variables);

        if (hasFine) {
          redistOps[lvl]->redistributeFine(*a_divG[lvl + 1], *m_massDifference[lvl], scale, variables);
        }
      }
    }
  }
}

void
CdrSolver::redistribute(EBAMRCellData& a_phi, const EBAMRIVData& a_delta) const noexcept
{
  CH_TIME("CdrSolver::redistribute");
  if (m_verbosity > 5) {
    pout() << m_name + "::redistribute" << endl;
  }

  CH_assert(a_phi.getRealm() == m_realm);
  CH_assert(a_delta.getRealm() == m_realm);

  Vector<RefCountedPtr<EBFluxRedistribution>>& redistOps = m_amr->getRedistributionOp(m_realm, m_phase);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Real     scale     = 1.0;
    const Interval variables = Interval(0, 0);
    const bool     hasCoar   = lvl > 0;
    const bool     hasFine   = lvl < m_amr->getFinestLevel();

    CH_assert(a_phi[lvl]->nComp() == 1);
    CH_assert(a_delta[lvl]->nComp() == 1);

    if (hasCoar) {
      redistOps[lvl]->redistributeCoar(*a_phi[lvl - 1], *a_delta[lvl], scale, variables);
    }

    redistOps[lvl]->redistributeLevel(*a_phi[lvl], *a_delta[lvl], scale, variables);

    if (hasFine) {
      redistOps[lvl]->redistributeFine(*a_phi[lvl + 1], *a_delta[lvl], scale, variables);
    }
  }
}

void
CdrSolver::computeAdvectionFlux(EBAMRFluxData&       a_flux,
                                const EBAMRFluxData& a_facePhi,
                                const EBAMRFluxData& a_faceVelocity,
                                const bool           a_addDomainFlux)
{
  CH_TIME("CdrSolver::computeAdvectionFlux(EBAMRFluxData, EBAMRFluxData, EBAMRFluxData, bool)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeAdvectionFlux(EBAMRFluxData, EBAMRFluxData, EBAMRFluxData, bool)" << endl;
  }

  CH_assert(a_flux[0]->nComp() == 1);
  CH_assert(a_facePhi[0]->nComp() == 1);
  CH_assert(a_faceVelocity[0]->nComp() == 1);

  // Computes the advection flux F = phi*v. This includes domain boundary faces.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    this->computeAdvectionFlux(*a_flux[lvl], *a_facePhi[lvl], *a_faceVelocity[lvl], lvl);
  }

  // Set domain faces to zero or enforce domain BCs (depending on flag)
  if (a_addDomainFlux) {
    this->fillDomainFlux(a_flux);
  }
  else {
    this->resetDomainFlux(a_flux);
  }
}

void
CdrSolver::computeAdvectionFlux(LevelData<EBFluxFAB>&       a_flux,
                                const LevelData<EBFluxFAB>& a_facePhi,
                                const LevelData<EBFluxFAB>& a_faceVelocity,
                                const int                   a_lvl)
{
  CH_TIME("CdrSolver::computeAdvectionFlux(LD<EBFluxFAB>, LD<EBFluxFAB>, LD<EBFluxFAB>, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeAdvectionFlux(LD<EBFluxFAB>, LD<EBFluxFAB>, LD<EBFluxFAB>, int)" << endl;
  }

  CH_assert(a_flux.nComp() == 1);
  CH_assert(a_facePhi.nComp() == 1);
  CH_assert(a_faceVelocity.nComp() == 1);

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const DataIterator&      dit = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    for (int dir = 0; dir < SpaceDim; dir++) {
      EBFaceFAB&       flux = a_flux[din][dir];
      const EBFaceFAB& phi  = a_facePhi[din][dir];
      const EBFaceFAB& vel  = a_faceVelocity[din][dir];

      flux.setVal(0.0, m_comp);
      flux += phi;
      flux *= vel;
    }
  }
}

void
CdrSolver::computeDiffusionFlux(EBAMRFluxData& a_flux, const EBAMRCellData& a_phi, const bool a_addDomainFlux)
{
  CH_TIME("CdrSolver::computeDiffusionFlux(EBAMRFluxData, EBAMRCellData, bool) ");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDiffusionFlux(EBAMRFluxData, EBAMRCellData, bool)" << endl;
  }

  CH_assert(a_flux[0]->nComp() == 1);
  CH_assert(a_phi[0]->nComp() == 1);

  // Computes the diffusion flux F = D*Grad(phi). This includes domain boundary faces.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    this->computeDiffusionFlux(*a_flux[lvl], *a_phi[lvl], lvl);
  }

  // Set domain faces to zero or enforce domain BCs (depending on flag)
  if (a_addDomainFlux) {
    this->fillDomainFlux(a_flux);
  }
  else {
    this->resetDomainFlux(a_flux);
  }
}

void
CdrSolver::computeDiffusionFlux(LevelData<EBFluxFAB>& a_flux, const LevelData<EBCellFAB>& a_phi, const int a_lvl)
{
  CH_TIME("CdrSolver::computeDiffusionFlux(LD<EBFluxFAB>, LD<EBCellFAB>, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDiffusionFlux(LD<EBFluxFAB>, LD<EBCellFAB>, int)" << endl;
  }

  CH_assert(a_flux.nComp() == 1);
  CH_assert(a_phi.nComp() == 1);

  // TLDR: This routine computes the diffusion flux F = D*Grad(phi) on face centers. Since this uses centered differencing
  //       and we don't have valid data outside the computational domain we only do the differencing for interior faces,
  //       setting the flux to zero on domain faces. Note that we need to fill flux in the tangential ghost face centers because
  //       the face centroid flux is interpolated between face centers.
  //

  const Real               dx        = m_amr->getDx()[a_lvl];
  const Real               inverseDx = 1. / dx;
  const DisjointBoxLayout& dbl       = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout&        ebisl     = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl];
  const ProblemDomain&     domain    = m_amr->getDomains()[a_lvl];
  const DataIterator&      dit       = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const EBCellFAB& phi     = a_phi[din];
    const Box&       cellBox = dbl[din];
    const EBISBox&   ebisbox = ebisl[din];
    const EBGraph&   ebgraph = ebisbox.getEBGraph();

    for (int dir = 0; dir < SpaceDim; dir++) {
      EBFaceFAB&       flux = a_flux[din][dir];
      const EBFaceFAB& dco  = (*m_faceCenteredDiffusionCoefficient[a_lvl])[din][dir];

      // Regular grid data.
      BaseFab<Real>&       regFlux = flux.getSingleValuedFAB();
      const BaseFab<Real>& regPhi  = phi.getSingleValuedFAB();
      const BaseFab<Real>& regDco  = dco.getSingleValuedFAB();

      flux.setVal(0.0);

      // Only want interior faces -- the domain flux will be set to zero (what else would it be...?). Anyways, recall that
      // the cut-cell face centroid fluxes are interpolated between face centers. So, if a face at the edge of a patch is cut
      // by the EB, the interpolant will reach out of the patch. Thus, we must fill the flux for the tangential ghost faces
      // outside the grid patch.
      Box grownCellBox = cellBox;
      grownCellBox.grow(1);
      grownCellBox &= domain;
      grownCellBox.grow(dir, -1);

      // These are the "regions" for the regular and cut-cell kernels.
      const Box    grownFaceBox = surroundingNodes(grownCellBox, dir);
      FaceIterator faceit(ebisbox.getIrregIVS(grownCellBox), ebgraph, dir, FaceStop::SurroundingWithBoundary);

      // Regular kernel. Note that we call the kernel on a face-centered box, so the cell on the high side is located at
      // iv, and the cell at the low side is at iv - BASISV(dir).
      auto regularKernel = [&](const IntVect& iv) -> void {
        regFlux(iv, m_comp) = inverseDx * regDco(iv, m_comp) * (regPhi(iv, m_comp) - regPhi(iv - BASISV(dir), m_comp));
      };

      // Cut-cell kernel. Basically the same as the above but we need to explicity get vofs on the low/high side (because
      // we may have multi-cells but the above kernel only does single-valued cells).
      auto irregularKernel = [&](const FaceIndex& face) -> void {
        if (!face.isBoundary()) {
          const VolIndex hiVoF = face.getVoF(Side::Hi);
          const VolIndex loVoF = face.getVoF(Side::Lo);

          flux(face, m_comp) = dco(face, m_comp) * (phi(hiVoF, m_comp) - phi(loVoF, m_comp)) / dx;
        }
      };

      // Execute kernels.
      BoxLoops::loop(grownFaceBox, regularKernel);
      BoxLoops::loop(faceit, irregularKernel);
    }
  }
}

void
CdrSolver::computeAdvectionDiffusionFlux(EBAMRFluxData&       a_flux,
                                         const EBAMRCellData& a_cellStates,
                                         const EBAMRFluxData& a_faceStates,
                                         const EBAMRFluxData& a_faceVelocities,
                                         const EBAMRFluxData& a_faceDiffCo,
                                         const bool           a_addDomainFlux)
{
  CH_TIME("CdrSolver::computeAdvectionDiffusionFlux(EBAMRFluxData, EBAMRCellData, EBAMRFluxDatax3, bool)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeAdvectionDiffusionFlux(EBAMRFluxData, EBAMRCellData, EBAMRFluxDatax3, bool)" << endl;
  }

  DataOps::setValue(a_flux, 0.0);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
    const ProblemDomain&     domain = m_amr->getDomains()[lvl];
    const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real               dx     = m_amr->getDx()[lvl];
    const Real               idx    = 1. / dx;
    const DataIterator&      dit    = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box&     cellBox = dbl[din];
      const EBISBox& ebisbox = ebisl[din];
      const EBGraph& ebgraph = ebisbox.getEBGraph();

      const EBCellFAB&     phiCell    = (*a_cellStates[lvl])[din];
      const BaseFab<Real>& regPhiCell = phiCell.getSingleValuedFAB();

      for (int dir = 0; dir < SpaceDim; dir++) {
        EBFaceFAB&       fluxFace = (*a_flux[lvl])[din][dir];
        const EBFaceFAB& phiFace  = (*a_faceStates[lvl])[din][dir];
        const EBFaceFAB& velFace  = (*a_faceVelocities[lvl])[din][dir];
        const EBFaceFAB& dcoFace  = (*a_faceDiffCo[lvl])[din][dir];

        // Regular grid data.
        BaseFab<Real>&       regFluxFace = fluxFace.getSingleValuedFAB();
        const BaseFab<Real>& regDcoFace  = dcoFace.getSingleValuedFAB();

        // Compute advective flux
        fluxFace.setVal(0.0, m_comp);
        fluxFace += phiFace;
        fluxFace *= velFace;

        const Box grownBox = grow(cellBox, 1) & domain;

        // When adding the diffusion flux we only do the interior faces (no diffusion flux across boundary faces). But
        // note that we need to fill fluxes also in the "ghost faces" because when computing the divergence we will
        // interpolate fluxes to centroids.
        Box interiorFaces = cellBox;
        interiorFaces.grow(1);
        interiorFaces &= domain;
        interiorFaces.grow(dir, -1);
        interiorFaces.surroundingNodes(dir);

        FaceIterator faceit(ebisbox.getIrregIVS(grownBox), ebgraph, dir, FaceStop::SurroundingNoBoundary);

        // Regular kernel. Note that we call the kernel on a face-centered box, so the cell on the high side is located at
        // iv, and the cell at the low side is at iv - BASISV(dir).
        auto interiorKernel = [&](const IntVect& iv) -> void {
          Real&       faceFlux  = regFluxFace(iv, m_comp);
          const Real& cellPhiLo = regPhiCell(iv - BASISV(dir), m_comp);
          const Real& cellPhiHi = regPhiCell(iv, m_comp);
          const Real& faceDco   = regDcoFace(iv, m_comp);

          faceFlux -= idx * faceDco * (cellPhiHi - cellPhiLo);
        };

        // Cut-cell kernel. Basically the same as the above but we need to explicity get vofs on the low/high side (because
        // we may have multi-cells but the above kernel only does single-valued cells).
        auto irregularKernel = [&](const FaceIndex& face) -> void {
          const VolIndex hiVoF = face.getVoF(Side::Hi);
          const VolIndex loVoF = face.getVoF(Side::Lo);

          const Real& cellPhiLo = phiCell(loVoF, m_comp);
          const Real& cellPhiHi = phiCell(hiVoF, m_comp);
          const Real& faceDco   = dcoFace(face, m_comp);

          fluxFace(face, m_comp) -= idx * faceDco * (cellPhiHi - cellPhiLo);
        };

        // Execute kernels.
        BoxLoops::loop(interiorFaces, interiorKernel);
        BoxLoops::loop(faceit, irregularKernel);
      }
    }
  }

  if (a_addDomainFlux) {
    this->fillDomainFlux(a_flux);
  }
}

void
CdrSolver::resetDomainFlux(EBAMRFluxData& a_flux)
{
  CH_TIME("CdrSolver::resetDomainFlux(EBAMRFluxData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::resetDomainFlux(EBAMRFluxData)" << endl;
  }

  CH_assert(a_flux[0]->nComp() == 1);

  // TLDR: This routine iterates through all faces in a_flux which are boundary faces and sets the flux there to zero.

  constexpr Real zero = 0.0;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
    const ProblemDomain&     domain = m_amr->getDomains()[lvl];
    const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const DataIterator&      dit    = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box      cellBox = dbl[din];
      const EBISBox& ebisbox = ebisl[din];
      const EBGraph& ebgraph = ebisbox.getEBGraph();

      for (int dir = 0; dir < SpaceDim; dir++) {
        EBFaceFAB&     flux    = (*a_flux[lvl])[din][dir];
        BaseFab<Real>& regFlux = flux.getSingleValuedFAB();

        for (SideIterator sit; sit.ok(); ++sit) {
          const Side::LoHiSide side = sit();

          // Create a cell box which lies next to the domain. Also make an iterator over the cut-cell faces.
          const Box    boundaryCellBox = adjCellBox(domain.domainBox(), dir, side, -1) & cellBox;
          FaceIterator faceit(ebisbox.getIrregIVS(boundaryCellBox), ebgraph, dir, FaceStop::AllBoundaryOnly);

          const int     sign  = (side == Side::Lo) ? 0 : 1;
          const IntVect shift = sign * BASISV(dir);

          // Regular kernel -- note the shift to make sure that cell indices map to face indices. I am doing this because
          // if we were to convert the box to a face-centered box we would do another layer of faces. So shift directly.
          auto regularKernel = [&](const IntVect& iv) -> void {
            regFlux(iv + shift, m_comp) = zero;
          };

          // Irregular kernel. Same as the above really.
          auto irregularKernel = [&](const FaceIndex& face) -> void {
            flux(face, m_comp) = 0.0;
          };

          // Execute kernels
          BoxLoops::loop(boundaryCellBox, regularKernel);
          BoxLoops::loop(faceit, irregularKernel);
        }
      }
    }
  }
}

void
CdrSolver::fillDomainFlux(EBAMRFluxData& a_flux)
{
  CH_TIME("CdrSolver::fillDomainFlux(EBAMRFluxData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::fillDomainFlux(EBAMRFluxData)" << endl;
  }

  CH_assert(a_flux[0]->nComp() == 1);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    this->fillDomainFlux(*a_flux[lvl], lvl);
  }
}

void
CdrSolver::fillDomainFlux(LevelData<EBFluxFAB>& a_flux, const int a_level)
{
  CH_TIME("CdrSolver::fillDomainFlux(LD<EBFluxFAB>, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::fillDomainFlux(LD<EBFluxFAB>, int)" << endl;
  }

  CH_assert(a_flux.nComp() == 1);

  // TLDR: This iterates through all domain faces and sets the BC flux to either data-based, function-based, wall bc, or extrapolated outflow. This
  //       routine uses a face-iterator (because of BaseIFFAB<T>), but performance should be acceptable since we go through a very small number of faces.

  const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[a_level];
  const ProblemDomain&     domain = m_amr->getDomains()[a_level];
  const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[a_level];
  const DataIterator&      dit    = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const EBISBox& ebisbox = ebisl[din];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const Box      cellBox = dbl[din];

    for (int dir = 0; dir < SpaceDim; dir++) {
      EBFaceFAB& flux = a_flux[din][dir];

      // Iterate over domain cells. I'm REALLY not sure about if this is performant...
      for (SideIterator sit; sit.ok(); ++sit) {
        const Side::LoHiSide side = sit();

        const CdrDomainBC::DomainSide    domainSide = std::make_pair(dir, side);
        const CdrDomainBC::BcType&       bcType     = m_domainBC.getBcType(domainSide);
        const CdrDomainBC::FluxFunction& bcFunction = m_domainBC.getBcFunction(domainSide);

        // Data-based BC holders -- needed if bcType is CdrDomainBC::BcType::DataBased.
        const BaseIFFAB<Real>& dataBasedFlux = (*m_domainFlux[a_level])[din](dir, side);

        // Create a box which abuts the current domain side
        Box boundaryCellBox;
        boundaryCellBox = adjCellBox(domain.domainBox(), dir, side, -1);
        boundaryCellBox &= cellBox;

        // Use a face-iterator here to go through domain faces -- sort of have to do this because of (the odd) design decision to use a BaseIFFAB for holding
        // the domain fluxes.
        //
        // A FaceIterator may be inefficient, but we are REALLY going through a very small number of faces so there shouldn't be a performance hit here.
        const IntVectSet ivs(boundaryCellBox);
        FaceIterator     faceit(ivs, ebgraph, dir, FaceStop::AllBoundaryOnly);

        auto kernel = [&](const FaceIndex& face) -> void {
          // Set the flux on the BC. This can be data-based, wall, function-based, outflow (extrapolated).
          switch (bcType) {
          case CdrDomainBC::BcType::DataBased: {
            flux(face, m_comp) = sign(sit()) * dataBasedFlux(face, m_comp);

            break;
          }
          case CdrDomainBC::BcType::Wall: {
            flux(face, m_comp) = 0.0;

            break;
          }
          case CdrDomainBC::BcType::Function: {
            const RealVect pos = EBArith::getFaceLocation(face, m_amr->getDx()[a_level], m_amr->getProbLo());

            flux(face, m_comp) = -sign(sit()) * bcFunction(pos, m_time);

            break;
          }
          case CdrDomainBC::BcType::Outflow: {
            flux(face, m_comp) = 0.0;

            // Get the next interior face(s) and extrapolate from these.
            const VolIndex interiorVof = face.getVoF(flip(sit()));

            const std::vector<FaceIndex> neighborFaces = ebisbox.getFaces(interiorVof, dir, flip(sit())).stdVector();

            if (neighborFaces.size() > 0) {
              Real sumArea = 0.0;

              for (const auto& f : neighborFaces) {
                const Real areaFrac = ebisbox.areaFrac(f);

                sumArea += areaFrac;
                flux(face, m_comp) += areaFrac * flux(f, m_comp);
              }

              flux(face, m_comp) = sign(sit()) * std::max((Real)0.0, sign(sit()) * flux(face, m_comp)) / sumArea;
            }

            break;
          }
          case CdrDomainBC::BcType::Solver: {
            // Don't do anything beacuse the solver will have filled the flux already.

            break;
          }
          default: {
            MayDay::Error(
              "CdrSolver::fillDomainFlux(LD<EBFluxFAB>, int) - trying to fill unsupported domain bc flux type");

            break;
          }
          }
        };
        // Run the kernel
        BoxLoops::loop(faceit, kernel);
      }
    }
  }
}

void
CdrSolver::conservativeDivergenceNoKappaDivision(EBAMRCellData&     a_conservativeDivergence,
                                                 EBAMRFluxData&     a_flux,
                                                 const EBAMRIVData& a_ebFlux)
{
  CH_TIME("CdrSolver::conservativeDivergenceNoKappaDivision(EBAMRCellData, EBAMRFluxData, EBAMRIVData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::conservativeDivergenceNoKappaDivision(EBAMRCellData, EBAMRFluxData, EBAMRIVData)" << endl;
  }

  CH_assert(a_conservativeDivergence[0]->nComp() == 1);
  CH_assert(a_flux[0]->nComp() == 1);
  CH_assert(a_ebFlux[0]->nComp() == 1);

  // This routine computes the "regular" finite volume conservative divergence Sum(fluxes) but does not divide by the volume fraction (which a naive implementation would).
  // Rather, this no-kappa-divided divergence is used for computing another divergence (a hybrid) divergence which represents a stable update to
  // the hyperbolic equation.
  //
  // This routine computes just that: a_conservativeDivergence = kappa*div(F) = Sum(fluxes).
  //
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    a_flux[lvl]->exchange();

    // Compute kappa*div(F) in regular cells
    this->conservativeDivergenceRegular(*a_conservativeDivergence[lvl], *a_flux[lvl], lvl);

    // Recompute divergence on irregular cells
    this->computeDivergenceIrregular(*a_conservativeDivergence[lvl], *a_flux[lvl], *a_ebFlux[lvl], lvl);

    a_conservativeDivergence[lvl]->exchange();
  }

  m_amr->conservativeAverage(a_conservativeDivergence, m_realm, m_phase);
  m_amr->interpGhost(a_conservativeDivergence, m_realm, m_phase);
}

void
CdrSolver::computeDivergenceIrregular(LevelData<EBCellFAB>&             a_divG,
                                      const LevelData<EBFluxFAB>&       a_centroidFlux,
                                      const LevelData<BaseIVFAB<Real>>& a_ebFlux,
                                      const int                         a_lvl)
{
  CH_TIME("CdrSolver::computeDivergenceIrregular(LD<EBCellFAB>, LD<EBFluxFAB>, LD<BaseIVFAB<Real> >, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDivergenceIrregular(LD<EBCellFAB>, LD<EBFluxFAB>, LD<BaseIVFAB<Real> >, int)" << endl;
  }

  CH_assert(a_divG.nComp() == 1);
  CH_assert(a_centroidFlux.nComp() == 1);
  CH_assert(a_ebFlux.nComp() == 1);

  // This computes the conservative divergence kappa*div(F) = sum(fluxes) in cut-cells. Fluxes
  // must be face centroid centered!

  const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl];
  const Real               dx    = m_amr->getDx()[a_lvl];
  const DataIterator&      dit   = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    EBCellFAB&             divG    = a_divG[din];
    const BaseIVFAB<Real>& ebflux  = a_ebFlux[din];
    const EBISBox&         ebisbox = ebisl[din];

    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[din];

    auto kernel = [&](const VolIndex& vof) -> void {
      const Real ebArea = ebisbox.bndryArea(vof);

      // Add the EB flux contribution to our sum(fluxes)
      divG(vof, m_comp) = ebflux(vof, m_comp) * ebArea;

      // Add face flux contributions to our sum(fluxes).
      for (int dir = 0; dir < SpaceDim; dir++) {
        const EBFaceFAB& flux = a_centroidFlux[din][dir];

        for (SideIterator sit; sit.ok(); ++sit) {
          const Side::LoHiSide side = sit();

          const int               isign = sign(side);
          const Vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, side);

          for (int iface = 0; iface < faces.size(); iface++) {
            const FaceIndex face     = faces[iface];
            const Real      faceArea = ebisbox.areaFrac(face);

            divG(vof, m_comp) += isign * faceArea * flux(face, m_comp);
          }
        }
      }

      // Scale divG by dx but not by kappa.
      divG(vof, m_comp) *= 1. / dx;
    };

    BoxLoops::loop(vofit, kernel);
  }
}

void
CdrSolver::conservativeDivergenceRegular(LevelData<EBCellFAB>&       a_divJ,
                                         const LevelData<EBFluxFAB>& a_flux,
                                         const int                   a_lvl)
{
  CH_TIME("CdrSolver::conservativeDivergenceRegular(LD<EBCellFAB>, LD<EBFluxFAB>, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::conservativeDivergenceRegular(LD<EBCellFAB>, LD<EBFluxFAB>, int)" << endl;
  }

  CH_assert(a_divJ.nComp() == 1);
  CH_assert(a_flux.nComp() == 1);

  // We compute the conservative divergence in regular grid cells. The divergence is in the form
  // kappa*div(J) = sum(fluxes)/dx

  const DisjointBoxLayout dbl       = m_amr->getGrids(m_realm)[a_lvl];
  const DataIterator&     dit       = dbl.dataIterator();
  const ProblemDomain     domain    = m_amr->getDomains()[a_lvl];
  const Real              dx        = m_amr->getDx()[a_lvl];
  const Real              inverseDx = 1. / dx;

  const int nbox = dit.size();
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const Box      cellBox = dbl[din];
    EBCellFAB&     divJ    = a_divJ[din];
    BaseFab<Real>& divJReg = divJ.getSingleValuedFAB();

    divJ.setVal(0.0);

    for (int dir = 0; dir < SpaceDim; dir++) {
      const EBFaceFAB&     flux    = a_flux[din][dir];
      const BaseFab<Real>& fluxReg = flux.getSingleValuedFAB();

      // Regular kernel. We call this for a cell-centered box so the high flux is on iv + BASISV(dir) and the low flux
      // on iv + BASISV(dir);
      auto regularKernel = [&](const IntVect& iv) -> void {
        divJReg(iv, m_comp) += inverseDx * (fluxReg(iv + BASISV(dir), m_comp) - fluxReg(iv, m_comp));
      };

      // Execute the kernel.
      BoxLoops::loop(cellBox, regularKernel);
    }

    // Reset irregular grid cells. These will be computed in a different way.
    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[din];

    auto irregularKernel = [&](const VolIndex& vof) -> void {
      divJ(vof, m_comp) = 0.0;
    };

    BoxLoops::loop(vofit, irregularKernel);
  }

  a_divJ.exchange();
}

void
CdrSolver::defineInterpolationStencils()
{
  CH_TIME("CdrSolver::defineInterpolationStencils()");
  if (m_verbosity > 5) {
    pout() << m_name + "::defineInterpolationStencils()" << endl;
  }

  // This routine defines interpolation stencils for going from face-centered fluxes to face-centroid fluxes.

  const int finestLevel = m_amr->getFinestLevel();

  for (int dir = 0; dir < SpaceDim; dir++) {
    (m_interpStencils[dir]).resize(1 + finestLevel);

    for (int lvl = 0; lvl <= finestLevel; lvl++) {
      const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
      const ProblemDomain&     domain = m_amr->getDomains()[lvl];
      const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const DataIterator&      dit    = dbl.dataIterator();

      m_interpStencils[dir][lvl] = RefCountedPtr<LayoutData<BaseIFFAB<FaceStencil>>>(
        new LayoutData<BaseIFFAB<FaceStencil>>(dbl));

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        BaseIFFAB<FaceStencil>& sten     = (*m_interpStencils[dir][lvl])[din];
        const Box               cellBox  = dbl[din];
        const EBISBox&          ebisbox  = ebisl[din];
        const EBGraph&          ebgraph  = ebisbox.getEBGraph();
        const IntVectSet        irregIVS = ebisbox.getIrregIVS(cellBox);

        sten.define(irregIVS, ebisbox.getEBGraph(), dir, m_nComp);

        FaceIterator faceIt(irregIVS, ebgraph, dir, FaceStop::SurroundingWithBoundary);

        auto kernel = [&](const FaceIndex& face) -> void {
          sten(face, m_comp) = EBArith::getInterpStencil(face, IntVectSet(), ebisbox, domain.domainBox());
        };

        BoxLoops::loop(faceIt, kernel);
      }
    }
  }
}

void
CdrSolver::initialData()
{
  CH_TIME("CdrSolver::initialData()");
  if (m_verbosity > 5) {
    pout() << m_name + "::initialData()" << endl;
  }

  // CdrSolver can be initialized in several way -- we can fill the initial data with analytic data from a mesh, or we can
  // deposit particles (with an NGP scheme) on the mesh. This function does both -- first particles if we have them and
  // then we increment with the function.

  DataOps::setValue(m_phi, 0.0);

  // Deposit particles if we have them.
  this->initialDataParticles();

  // Increment with function values if this is also called for. Note that this increments,
  // why is why initialDataParticles is called first!
  this->initialDataDistribution();

  m_amr->conservativeAverage(m_phi, m_realm, m_phase);
  m_amr->interpGhost(m_phi, m_realm, m_phase);
}

void
CdrSolver::initialDataDistribution()
{
  CH_TIME("CdrSolver::initialDataDistribution()");
  if (m_verbosity > 5) {
    pout() << m_name + "::initialDataDistribution()" << endl;
  }

  // TLDR: We just run through every cell in the grid and increment by m_species->initialData.

  // Expose the initial data function to something DataOps can use.
  auto initFunc = [&](const RealVect& pos) -> Real {
    return m_species->initialData(pos, m_time);
  };

  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, m_nComp);

  DataOps::setValue(scratch, initFunc, m_amr->getProbLo(), m_amr->getDx(), m_comp);
  DataOps::incr(m_phi, scratch, 1.0);
  DataOps::setCoveredValue(m_phi, 0, 0.0);

  m_amr->conservativeAverage(m_phi, m_realm, m_phase);
  m_amr->interpGhost(m_phi, m_realm, m_phase);
}

void
CdrSolver::initialDataParticles()
{
  CH_TIME("CdrSolver::initialDataParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::initialDataParticles" << endl;
  }

  // TLDR: This function deposits a list of initial particles on the mesh using
  //       an NGP scheme, ignoring conservation in cut-cells. Please refers to
  //       the particle API to see the details of the ParticleContainer<T> type and
  //       deposition methods.

  const List<PointParticle>& initialParticles = m_species->getInitialParticles();

  const long long numParticles = (long long)initialParticles.length();

  if (ParallelOps::sum(numParticles) > 0LL) {

    // Make a ParticleContainer<T> and redistribute particles over the AMR hierarchy.
    ParticleContainer<PointParticle> particles;
    m_amr->allocate(particles, m_realm);
    particles.addParticles(m_species->getInitialParticles());

    // This function will be called BEFORE initialDataFunction, we it is safe to set m_phi to zero
    DataOps::setValue(m_phi, 0.0);

    // Deposit onto mesh.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const RealVect           dx     = m_amr->getDx()[lvl] * RealVect::Unit;
      const RealVect           probLo = m_amr->getProbLo();
      const ProblemDomain&     domain = m_amr->getDomains()[lvl];
      const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const DataIterator&      dit    = dbl.dataIterator();

      // 2. Deposit this levels particles. We use an NGP scheme to do this.
      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        const Box      cellBox = dbl[din];
        const EBISBox& ebisbox = ebisl[din];

        // Make the deposition object and put the particles on the grid.
        constexpr bool forceIrregNGP = true;
        EBParticleMesh interp(domain, cellBox, ebisbox, dx, probLo);

        interp.deposit<PointParticle, &PointParticle::weight>(particles[lvl][din].listItems(),
                                                              (*m_phi[lvl])[din],
                                                              DepositionType::NGP,
                                                              forceIrregNGP);
      }
    }
  }
}

void
CdrSolver::hybridDivergence(EBAMRCellData&     a_hybridDivergence,
                            EBAMRIVData&       a_massDifference,
                            const EBAMRIVData& a_nonConservativeDivergence)
{
  CH_TIME("CdrSolver::hybridDivergence(EBAMRCellData, EBAMRIVData, EBAMRIVData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::hybridDivergence(EBAMRCellData, EBAMRIVData, EBAMRIVData)" << endl;
  }

  CH_assert(a_hybridDivergence[0]->nComp() == 1);
  CH_assert(a_massDifference[0]->nComp() == 1);
  CH_assert(a_nonConservativeDivergence[0]->nComp() == 1);

  // On input, a_hybridDivergence must contain kappa*div(F). We also have the non-conservative divergence computed, so we just compute
  // the regular hybrid divergence and the mass loss/gain.

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    this->hybridDivergence(*a_hybridDivergence[lvl], *a_massDifference[lvl], *a_nonConservativeDivergence[lvl], lvl);
  }
}

void
CdrSolver::hybridDivergence(LevelData<EBCellFAB>&             a_hybridDivergence,
                            LevelData<BaseIVFAB<Real>>&       a_massDifference,
                            const LevelData<BaseIVFAB<Real>>& a_nonConservativeDivergence,
                            const int                         a_lvl)
{
  CH_TIME("CdrSolver::hybridDivergence(LD<EBCellFAB>, LD<BaseIVFAB<Real> >, LD<BaseIVFAB<Real> >, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::hybridDivergence(LD<EBCellFAB>, LD<BaseIVFAB<Real> >, LD<BaseIVFAB<Real> >, int)" << endl;
  }

  CH_assert(a_hybridDivergence.nComp() == 1);
  CH_assert(a_massDifference.nComp() == 1);
  CH_assert(a_nonConservativeDivergence.nComp() == 1);

  // TLDR: We want to compute a stable approximation to the divergence using Dh = kappa*Dc + (1-kappa)*Dnc where
  //       Dh is the hybrid divergence and Dc/Dnc are the conservative and non-conservative divergences. On the way in
  //       we had a_hybridDivergence = kappa*Dc.

  const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl];
  const DataIterator&      dit   = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    // On input, divH contains kappa*div(F) and divNC contains the non-conservative divergence
    EBCellFAB&             divH   = a_hybridDivergence[din];
    BaseIVFAB<Real>&       deltaM = a_massDifference[din];
    const BaseIVFAB<Real>& divNC  = a_nonConservativeDivergence[din];

    const EBISBox& ebisbox = ebisl[din];

    VoFIterator& vofit  = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[din];
    auto         kernel = [&](const VolIndex& vof) -> void {
      const Real kappa = ebisbox.volFrac(vof);
      const Real dc    = divH(vof, m_comp);
      const Real dnc   = divNC(vof, m_comp);

      // Note to self: deltaM = (1-kappa)*(dc - kappa*dnc) because dc was not divided by kappa,
      // which it would be otherwise.
      divH(vof, m_comp)   = dc + (1 - kappa) * dnc;           // On output, contains hybrid divergence
      deltaM(vof, m_comp) = (1 - kappa) * (dc - kappa * dnc); // On output, contains mass loss/gain.
    };

    BoxLoops::loop(vofit, kernel);
  }
}

void
CdrSolver::interpolateFluxToFaceCentroids(EBAMRFluxData& a_flux)
{
  CH_TIME("CdrSolver::interpolateFluxToFaceCentroids(EBAMRFluxData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::interpolateFluxToFaceCentroids(EBAMRFluxData)" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    a_flux[lvl]->exchange();

    this->interpolateFluxToFaceCentroids(*a_flux[lvl], lvl);
    ;
  }
}

void
CdrSolver::interpolateFluxToFaceCentroids(LevelData<EBFluxFAB>& a_flux, const int a_lvl)
{
  CH_TIME("CdrSolver::interpolateFluxToFaceCentroids(LD<EBFluxFAB>, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::interpolateFluxToFaceCentroids(LD<EBFluxFAB>, int)" << endl;
  }

  CH_assert(a_flux.nComp() == 1);

  // TLDR: We are given face-centered fluxes which we want to put on face centroids. To do this we store face-centered
  //       fluxes on temporary storage and just interpolate them to face centroids.

  const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl];
  const DataIterator&      dit   = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const Box        cellBox  = dbl.get(din);
    const EBISBox&   ebisbox  = ebisl[din];
    const EBGraph&   ebgraph  = ebisbox.getEBGraph();
    const IntVectSet irregIVS = ebisbox.getIrregIVS(cellBox);

    const bool isRegular   = ebisbox.isRegular(cellBox);
    const bool isCovered   = ebisbox.isCovered(cellBox);
    const bool isIrregular = !isRegular && !isCovered;

    if (isIrregular) {
      for (int dir = 0; dir < SpaceDim; dir++) {
        EBFaceFAB& faceFlux = a_flux[din][dir];

        // Make a copy of faceFlux which we use as interpolant. EBFaceFAB
        // needs a cell box, so here you go.
        EBFaceFAB clone(ebisbox, faceFlux.getCellRegion(), dir, m_nComp);
        clone.setVal(0.0);
        clone += faceFlux;

        // Compute face centroid flux on cut-cell face centroids. Since a_flux enforces boundary conditions
        // we include domain boundary cut-cell faces in the interpolation.
        FaceIterator faceit(irregIVS, ebgraph, dir, FaceStop::SurroundingWithBoundary);

        auto kernel = [&](const FaceIndex& face) -> void {
          const FaceStencil& sten = (*m_interpStencils[dir][a_lvl])[din](face, m_comp);

          faceFlux(face, m_comp) = 0.;

          for (int i = 0; i < sten.size(); i++) {
            const FaceIndex& iface   = sten.face(i);
            const Real       iweight = sten.weight(i);

            faceFlux(face, m_comp) += iweight * clone(iface, m_comp);
          }
        };

        BoxLoops::loop(faceit, kernel);
      }
    }
  }
}

void
CdrSolver::nonConservativeDivergence(EBAMRIVData& a_nonConservativeDivergence, const EBAMRCellData& a_divG)
{
  CH_TIME("CdrSolver::nonConservativeDivergence(EBAMRIVData, EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::nonConservativeDivergence(EBAMRIVData, EBAMRCellData)" << endl;
  }

  CH_assert(a_nonConservativeDivergence[0]->nComp() == 1);
  CH_assert(a_divG[0]->nComp() == 1);

  // TLDR: This routine computes the non-conservative divergence divNC(G) = sum(kappa*div(G))/sum(kappa). This is done through
  //       AmrMesh's stencil in cut cells. The neighborhood consists of cells that can be reached with a monotone path.

  if (m_blendConservation) {
    m_amr->nonConservativeDivergence(a_nonConservativeDivergence, a_divG, m_realm, m_phase);
  }
  else { // Mostly for debugging.
    DataOps::setValue(a_nonConservativeDivergence, 0.0);
  }
}

void
CdrSolver::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("CdrSolver::regrid(int, int, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::regrid(int, int, int)" << endl;
  }

  // Allocate internal storage. This will fill m_phi with invalid data, but preRegrid(...) should have been called
  // prior to this routine so the old m_phi (before the regrid) is stored in m_cachedPhi.
  this->allocate();

  const EBCoarseToFineInterp::Type interpType = m_regridSlopes ? EBCoarseToFineInterp::Type::ConservativeMinMod
                                                               : EBCoarseToFineInterp::Type::ConservativePWC;

  // Interpolate to the new grids.
  m_amr->interpToNewGrids(m_phi, m_cachePhi, m_phase, a_lmin, a_oldFinestLevel, a_newFinestLevel, interpType);
  m_amr->interpToNewGrids(m_source, m_cacheSource, m_phase, a_lmin, a_oldFinestLevel, a_newFinestLevel, interpType);

  // Coarsen data and update ghost cells.
  m_amr->conservativeAverage(m_phi, m_realm, m_phase);
  m_amr->interpGhost(m_phi, m_realm, m_phase);

  m_amr->conservativeAverage(m_source, m_realm, m_phase);
  m_amr->interpGhost(m_source, m_realm, m_phase);

  // Deallocate the scratch storage.
  m_cachePhi.clear();
  m_cacheSource.clear();
}

void
CdrSolver::setAmr(const RefCountedPtr<AmrMesh>& a_amr)
{
  CH_TIME("CdrSolver::setAmr(AmrMesh)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setAmr(AmrMesh)" << endl;
  }

  CH_assert(!a_amr.isNull());

  m_amr = a_amr;
}

void
CdrSolver::registerOperators()
{
  CH_TIME("CdrSolver::registerOperators()");
  if (m_verbosity > 5) {
    pout() << m_name + "::registerOperators()" << endl;
  }

  CH_assert(!m_amr.isNull());

  m_amr->registerOperator(s_eb_coar_ave, m_realm, m_phase);
  m_amr->registerOperator(s_eb_fill_patch, m_realm, m_phase);
  m_amr->registerOperator(s_eb_fine_interp, m_realm, m_phase);
  m_amr->registerOperator(s_eb_redist, m_realm, m_phase);
  m_amr->registerOperator(s_eb_irreg_interp, m_realm, m_phase);
  m_amr->registerOperator(s_noncons_div, m_realm, m_phase);
}

void
CdrSolver::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry)
{
  CH_TIME("CdrSolver::setComputationalGeometry(RefCountedPtr<ComputationalGeometry>)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setComputationalGeometry(RefCountedPtr<ComputationalGeometry>)" << endl;
  }

  CH_assert(!a_computationalGeometry.isNull());

  m_computationalGeometry = a_computationalGeometry;

  const RefCountedPtr<MultiFluidIndexSpace> MultiFluidIndexSpace = m_computationalGeometry->getMfIndexSpace();
}

void
CdrSolver::setDiffusionCoefficient(const EBAMRFluxData& a_diffusionCoefficient,
                                   const EBAMRIVData&   a_ebDiffusionCoefficient)
{
  CH_TIME("CdrSolver::setDiffusionCoefficient(EBAMRFluxData, EBAMRIVData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setDiffusionCoefficient(EBAMRFluxData, EBAMRIVData)" << endl;
  }

  CH_assert(a_diffusionCoefficient[0]->nComp() == 1);
  CH_assert(a_ebDiffusionCoefficient[0]->nComp() == 1);

  // Do a copy -- realms do not have to be the same.
  m_amr->copyData(m_faceCenteredDiffusionCoefficient, a_diffusionCoefficient);
  m_amr->copyData(m_ebCenteredDiffusionCoefficient, a_ebDiffusionCoefficient);

  m_amr->conservativeAverage(m_faceCenteredDiffusionCoefficient, m_realm, m_phase);
  m_amr->conservativeAverage(m_ebCenteredDiffusionCoefficient, m_realm, m_phase);
}

void
CdrSolver::setDiffusionCoefficient(const Real a_diffusionCoefficient)
{
  CH_TIME("CdrSolver::setDiffusionCoefficient(Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setDiffusionCoefficient(Real)" << endl;
  }

  DataOps::setValue(m_cellCenteredDiffusionCoefficient, a_diffusionCoefficient);
  DataOps::setValue(m_faceCenteredDiffusionCoefficient, a_diffusionCoefficient);
  DataOps::setValue(m_ebCenteredDiffusionCoefficient, a_diffusionCoefficient);

  m_amr->conservativeAverage(m_cellCenteredDiffusionCoefficient, m_realm, m_phase);
  m_amr->conservativeAverage(m_faceCenteredDiffusionCoefficient, m_realm, m_phase);
  m_amr->conservativeAverage(m_ebCenteredDiffusionCoefficient, m_realm, m_phase);
}

void
CdrSolver::setDiffusionCoefficient(const std::function<Real(const RealVect a_position)>& a_diffCo)
{
  CH_TIME("CdrSolver::setDiffusionCoefficient(std::function<Real(const RealVect a_position)>))");
  if (m_verbosity > 5) {
    pout() << m_name + "::setDiffusionCoefficient(std::function<Real(const RealVect a_position)>)" << endl;
  }

  DataOps::setValue(m_cellCenteredDiffusionCoefficient, a_diffCo, m_amr->getProbLo(), m_amr->getDx(), m_comp);
  DataOps::setValue(m_faceCenteredDiffusionCoefficient, a_diffCo, m_amr->getProbLo(), m_amr->getDx(), m_comp);
  DataOps::setValue(m_ebCenteredDiffusionCoefficient, a_diffCo, m_amr->getProbLo(), m_amr->getDx(), m_comp);
}

void
CdrSolver::setEbFlux(const EBAMRIVData& a_ebFlux)
{
  CH_TIME("CdrSolver::setEbFlux(a_ebFlux)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setEbFlux(a_ebFlux)" << endl;
  }

  CH_assert(a_ebFlux[0]->nComp() == 1);

  m_amr->copyData(m_ebFlux, a_ebFlux);
}

void
CdrSolver::setEbFlux(const Real a_ebFlux)
{
  CH_TIME("CdrSolver::setEbFlux(Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setEbFlux(Real)" << endl;
  }

  DataOps::setValue(m_ebFlux, a_ebFlux);
}

void
CdrSolver::setSpecies(const RefCountedPtr<CdrSpecies>& a_species)
{
  CH_TIME("CdrSolver::setSpecies(RefCountedPtr<CdrSpecies>)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setSpecies(RefCountedPtr<CdrSpecies>)" << endl;
  }

  CH_assert(!a_species.isNull());

  m_species     = a_species;
  m_name        = m_species->getName();
  m_isDiffusive = m_species->isDiffusive();
  m_isMobile    = m_species->isMobile();
}

void
CdrSolver::setSource(const EBAMRCellData& a_source)
{
  CH_TIME("CdrSolver::setSource(EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setSource(EBAMRCellData)" << endl;
  }

  CH_assert(a_source[0]->nComp() == 1);

  m_amr->copyData(m_source, a_source);

  m_amr->conservativeAverage(m_source, m_realm, m_phase);
  m_amr->interpGhost(m_source, m_realm, m_phase);
}

void
CdrSolver::setSource(const Real a_source)
{
  CH_TIME("CdrSolver::setSource(Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setSource(Real)" << endl;
  }

  DataOps::setValue(m_source, a_source);

  m_amr->conservativeAverage(m_source, m_realm, m_phase);
  m_amr->interpGhost(m_source, m_realm, m_phase);
}

void
CdrSolver::setSource(const std::function<Real(const RealVect a_position)> a_source)
{
  CH_TIME("CdrSolver::setSource(std::function<Real(const RealVect a_position)>)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setSource(std::function<Real(const RealVect a_position)>)" << endl;
  }

  DataOps::setValue(m_source, a_source, m_amr->getProbLo(), m_amr->getDx(), m_comp);

  m_amr->conservativeAverage(m_source, m_realm, m_phase);
  m_amr->interpGhost(m_source, m_realm, m_phase);
}

void
CdrSolver::setTime(const int a_step, const Real a_time, const Real a_dt)
{
  CH_TIME("CdrSolver::setTime(int, Real, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setTime(int, Real, Real)" << endl;
  }

  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;
}

void
CdrSolver::setVelocity(const EBAMRCellData& a_velo)
{
  CH_TIME("CdrSolver::setVelocity");
  if (m_verbosity > 5) {
    pout() << m_name + "::setVelocity" << endl;
  }

  CH_assert(a_velo[0]->nComp() == SpaceDim);

  m_amr->copyData(m_cellVelocity, a_velo);

  m_amr->conservativeAverage(m_cellVelocity, m_realm, m_phase);
  m_amr->interpGhost(m_cellVelocity, m_realm, m_phase);
}

void
CdrSolver::setVelocity(const RealVect a_velo)
{
  CH_TIME("CdrSolver::setVelocity(RealVect)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setVelocity(RealVect)" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    for (int dir = 0; dir < SpaceDim; dir++) {
      DataOps::setValue(*m_cellVelocity[lvl], a_velo[dir], dir);
    }
  }

  m_amr->conservativeAverage(m_cellVelocity, m_realm, m_phase);
  m_amr->interpGhost(m_cellVelocity, m_realm, m_phase);
}

void
CdrSolver::setVelocity(const std::function<RealVect(const RealVect a_pos)>& a_velo)
{
  CH_TIME("CdrSolver::setVelocity(std::function<RealVect(const RealVect a_pos)>)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setVelocity(std::function<RealVect(const RealVect a_pos)>)" << endl;
  }

  DataOps::setValue(m_cellVelocity, a_velo, m_amr->getProbLo(), m_amr->getDx());
}

void
CdrSolver::setPhase(const phase::which_phase a_phase)
{
  CH_TIME("CdrSolver::setPhase(phase::which_phase)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setPhase(phase::which_phase)" << endl;
  }

  m_phase = a_phase;
}

void
CdrSolver::setVerbosity(const int a_verbosity)
{
  CH_TIME("CdrSolver::setVerbosity(int)");
  m_verbosity = a_verbosity;

  if (m_verbosity > 5) {
    pout() << m_name + "::setVerbosity(int)" << endl;
  }
}

void
CdrSolver::writePlotFile()
{
  CH_TIME("CdrSolver::writePlotFile()");
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotFile()" << endl;
  }

  // We recycle the writePlotData(...) routine to write a plot file. That routine is used for external
  // calls to this solver when we want to populate an external plot file with data from this class. Nothing
  // wrong with reusing it here.

  // Number of output components and their names
  const int                 numPlotVars  = this->getNumberOfPlotVariables();
  const Vector<std::string> plotVarNames = this->getPlotVariableNames();

  // Allocate storage
  EBAMRCellData output;
  m_amr->allocate(output, m_realm, m_phase, numPlotVars);
  DataOps::setValue(output, 0.0);

  // Copy internal data to be plotted over to 'output'
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    int icomp = 0;
    this->writePlotData(*output[lvl], icomp, m_realm, lvl);
  }

  // Filename
  char filename[100];
  sprintf(filename, "%s.step%07d.%dd.hdf5", m_name.c_str(), m_timeStep, SpaceDim);
  std::string fname(filename);

  // Alias, because Chombo's EBAMRIO wants raw pointers (but we stick to smart pointers, like God intended).
  Vector<LevelData<EBCellFAB>*> outputPtr;
  m_amr->alias(outputPtr, output);

#ifdef CH_USE_HDF5
  constexpr int numPlotGhost = 0;

  DischargeIO::writeEBHDF5(fname,
                           plotVarNames,
                           m_amr->getGrids(m_realm),
                           outputPtr,
                           m_amr->getDomains(),
                           m_amr->getDx(),
                           m_amr->getRefinementRatios(),
                           m_dt,
                           m_time,
                           m_amr->getProbLo(),
                           m_amr->getFinestLevel() + 1,
                           numPlotGhost);
#endif
}

void
CdrSolver::writePlotData(LevelData<EBCellFAB>& a_output,
                         int&                  a_icomp,
                         const std::string     a_outputRealm,
                         const int             a_level) const noexcept
{
  CH_TIME("CdrSolver::writePlotData");
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotData" << endl;
  }

  // TLDR: This routine writes plot data to an "external" data holder a_output, starting on a_icomp. The intention
  //       behind that is that a_output is a plot file which accumulates data over many solvers. Note that because
  //       this solver is cell-centered, we interpolate m_phi to centroids when outputting it.

  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, 1);

  // Plot state
  if (m_plotPhi) {
    this->writeData(a_output, a_icomp, m_phi, a_outputRealm, a_level, true, true);
    if (m_plotNumbers) {
      LevelData<EBCellFAB> alias;

      aliasLevelData<EBCellFAB>(alias, &a_output, Interval(a_icomp - 1, a_icomp - 1));

      DataOps::scale(alias, std::pow(m_amr->getDx()[a_level], SpaceDim));
    }
  }

  // Plot diffusion coefficients. These are stored on face centers but we need cell-centered data for output.
  if (m_plotDiffusionCoefficient && m_isDiffusive) {
    DataOps::averageFaceToCell(*scratch[a_level],
                               *m_faceCenteredDiffusionCoefficient[a_level],
                               m_amr->getDomains()[a_level]);

    // Do the previous because we need the ghost cells too.
    if (a_level > 0) {
      DataOps::averageFaceToCell(*scratch[a_level - 1],
                                 *m_faceCenteredDiffusionCoefficient[a_level - 1],
                                 m_amr->getDomains()[a_level - 1]);
    }

    this->writeData(a_output, a_icomp, scratch, a_outputRealm, a_level, false, true);
  }

  // Plot source terms
  if (m_plotSource) {
    this->writeData(a_output, a_icomp, m_source, a_outputRealm, a_level, false, true);

    if (m_plotNumbers) {
      LevelData<EBCellFAB> alias;

      aliasLevelData<EBCellFAB>(alias, &a_output, Interval(a_icomp - 1, a_icomp - 1));

      DataOps::scale(alias, std::pow(m_amr->getDx()[a_level], SpaceDim));
    }
  }

  // Plot velocities
  if (m_plotVelocity && m_isMobile) {
    this->writeData(a_output, a_icomp, m_cellVelocity, a_outputRealm, a_level, false, true);
  }

  // Plot EB fluxes. These are stored on sparse data structures but we need them on cell centers. So copy them to scratch and write data.
  if (m_plotEbFlux && m_isMobile) {
    DataOps::setValue(*scratch[a_level], 0.0);
    DataOps::incr(*scratch[a_level], *m_ebFlux[a_level], 1.0);
    this->writeData(a_output, a_icomp, scratch, a_outputRealm, a_level, false, false);
  }
}

void
CdrSolver::writeData(LevelData<EBCellFAB>& a_output,
                     int&                  a_comp,
                     const EBAMRCellData&  a_data,
                     const std::string     a_outputRealm,
                     const int             a_level,
                     const bool            a_interpToCentroids,
                     const bool            a_interpGhost) const noexcept

{
  CH_TIMERS("CdrSolver::writeData");
  CH_TIMER("CdrSolver::writeData::allocate", t1);
  CH_TIMER("CdrSolver::writeData::local_copy", t2);
  CH_TIMER("CdrSolver::writeData::interp_ghost", t3);
  CH_TIMER("CdrSolver::writeData::interp_centroid", t4);
  CH_TIMER("CdrSolver::writeData::final_copy", t5);
  if (m_verbosity > 5) {
    pout() << m_name + "::writeData" << endl;
  }

  // Number of components we are working with.
  const int numComp = a_data[a_level]->nComp();

  // Component ranges that we copy to/from.
  const Interval srcInterv(0, numComp - 1);
  const Interval dstInterv(a_comp, a_comp + numComp - 1);

  CH_START(t1);
  LevelData<EBCellFAB> scratch;
  m_amr->allocate(scratch, m_realm, m_phase, a_level, numComp);
  CH_STOP(t1);

  CH_START(t2);
  m_amr->copyData(scratch, *a_data[a_level], a_level, m_realm, m_realm);
  CH_START(t2);

  // Interpolate ghost cells
  CH_START(t3);
  if (a_level > 0 && a_interpGhost) {
    m_amr->interpGhost(scratch, *a_data[a_level - 1], a_level, m_realm, m_phase);
  }
  CH_STOP(t3);

  CH_START(t4);
  if (a_interpToCentroids) {
    m_amr->interpToCentroids(scratch, m_realm, m_phase, a_level);
  }
  CH_STOP(t4);

  DataOps::setCoveredValue(scratch, 0.0);

  CH_START(t5);
  m_amr->copyData(a_output, scratch, a_level, a_outputRealm, m_realm, dstInterv, srcInterv);
  CH_STOP(t5);

  a_comp += numComp;
}

#ifdef CH_USE_HDF5
void
CdrSolver::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const
{
  CH_TIME("CdrSolver::writeCheckpointLevel(HDF5Handle, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::writeCheckpointLevel(HDF5Handle, int)" << endl;
  }

  // Write state vector
  write(a_handle, *m_phi[a_level], m_name);
  write(a_handle, *m_source[a_level], m_name + "_src");
}
#endif

#ifdef CH_USE_HDF5
void
CdrSolver::readCheckpointLevel(HDF5Handle& a_handle, const int a_level)
{
  CH_TIME("CdrSolver::readCheckpointLevel(HDF5Handle, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::readCheckpointLevel(HDF5Handle, int)" << endl;
  }

  const Interval interv(m_comp, m_comp);

  read<EBCellFAB>(a_handle, *m_phi[a_level], m_name, m_amr->getGrids(m_realm)[a_level], interv, false);
  read<EBCellFAB>(a_handle, *m_source[a_level], m_name + "_src", m_amr->getGrids(m_realm)[a_level], interv, false);
}
#endif

Real
CdrSolver::computeAdvectionDt()
{
  CH_TIME("CdrSolver::computeAdvectionDt()");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeAdvectionDt()" << endl;
  }

  // TLDR: For advection we must have dt <= dx/(|vx|+|vy|+|vz|). E.g., with first order upwind phi^(k+1)_i = phi^k_i - (v*dt) * (phi^k_i - phi^k_(i-1))/dx so
  //       if phi^k_(i-1) == 0 then (1 - v*dt/dx) > 0.0 yields a positive definite solution (more general analysis when we have limiters is probably possible...)

  Real minDt = std::numeric_limits<Real>::max();

  if (m_isMobile) {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const Real               dx    = m_amr->getDx()[lvl];
      const DataIterator&      dit   = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(min : minDt)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        const Box        cellBox = dbl[din];
        const EBCellFAB& velo    = (*m_cellVelocity[lvl])[din];
        const EBISBox&   ebisBox = ebisl[din];

        VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[din];

        // Regular grid data.
        const BaseFab<Real>& veloReg = velo.getSingleValuedFAB();

        // Compute dt = dx/(|vx|+|vy|+|vz|) and check if it's smaller than the smallest so far.
        auto regularKernel = [&](const IntVect& iv) -> void {
          if (!ebisBox.isCovered(iv)) {
            Real vel = 0.0;
            for (int dir = 0; dir < SpaceDim; dir++) {
              vel += std::abs(veloReg(iv, dir));
            }

            minDt = std::min(dx / vel, minDt);
          }
        };

        // Same kernel, but for cut-cells.
        auto irregularKernel = [&](const VolIndex& vof) -> void {
          Real vel = 0.0;
          for (int dir = 0; dir < SpaceDim; dir++) {
            vel += std::abs(velo(vof, dir));
          }

          minDt = std::min(dx / vel, minDt);
        };

        // Execute the kernels.
        BoxLoops::loop(cellBox, regularKernel);
        BoxLoops::loop(vofit, irregularKernel);
      }
    }
  }

  return ParallelOps::min(minDt);
}

Real
CdrSolver::computeDiffusionDt()
{
  CH_TIME("CdrSolver::computeDiffusionDt()");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDiffusionDt()" << endl;
  }

  // TLDR: For advection we must have dt <= (dx*dx)/(2*d*D) where D is diffusion coefficient and d is spatial dimensions.
  //       E.g. in 1D, centered differencing yields
  //
  //          phi^(k+1)_i = phi^k_i - dt*D*(phi^k_(i+1) - 2*phi^k_i + phi^k_(i-1))/(dx*dx).
  //
  //       So for phi^k_(i-1) = phi^k_(i+1) = 0 we have phi^(k+1)_i = phi^k * (1 - 2*dt*D/(dx*dx)) which is positive definite only for dt < (dx*dx)/(2*D). Again,
  //       a more general analysis could be made.

  Real minDt = std::numeric_limits<Real>::max();

  if (m_isDiffusive) {

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const Real               dx    = m_amr->getDx()[lvl];
      const Real               dx2   = dx * dx;
      const DataIterator&      dit   = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(min : minDt)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        const Box        cellBox    = dbl[din];
        const EBISBox&   ebisbox    = ebisl[din];
        const EBFluxFAB& diffCoFace = (*m_faceCenteredDiffusionCoefficient[lvl])[din];
        VoFIterator&     vofit      = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[din];

        // Regular kernel. Strictly speaking, we should have a kernel which increments with the diffusion coefficients on each face since the finite volume
        // approximation to the Laplacian becomes (in 1D) Div*(D*Grad(phi)) = -D_(i-1/2)*(phi_i - phi_(i-1)) + D_(i+1/2)*(phi_(i+1)-phi_i). A good kernel
        // would do just that, but a lazy programmer just find the largest diffusion coefficient and uses that as an approximation.
        for (int dir = 0; dir < SpaceDim; dir++) {
          const BaseFab<Real>& diffCoReg = diffCoFace[dir].getSingleValuedFAB();

          auto regularKernel = [&](const IntVect& iv) -> void {
            if (ebisbox.isRegular(iv)) {

              const Real loD = diffCoReg(iv, m_comp);
              const Real hiD = diffCoReg(iv + BASISV(dir), m_comp);

              minDt = std::min(minDt, dx2 / (2 * SpaceDim * std::max(loD, hiD)));
            }
          };

          // Execute the kernel.
          BoxLoops::loop(cellBox, regularKernel);
        }

        // Same kernel as above, but we need to fetch grid faces differently.
        auto irregularKernel = [&](const VolIndex& vof) -> void {
          Real D = std::numeric_limits<Real>::min();

          for (int dir = 0; dir < SpaceDim; dir++) {
            const EBFaceFAB& Dface = diffCoFace[dir];

            for (SideIterator sit; sit.ok(); ++sit) {
              const std::vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, sit()).stdVector();

              for (const auto& face : faces) {
                D = std::max(D, Dface(face, m_comp));
              }
            }
          }

          minDt = std::min(minDt, dx2 / (2 * SpaceDim * D));
        };

        // Execute the kernel
        BoxLoops::loop(vofit, irregularKernel);
      }
    }
  }

  return ParallelOps::min(minDt);
}

Real
CdrSolver::computeAdvectionDiffusionDt()
{
  CH_TIME("CdrSolver::computeAdvectionDiffusionDt()");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeAdvectionDiffusionDt()" << endl;
  }

  // In 1D we have, e.g. d(phi)/dt = -d/dx(v*phi) + D*d^2(phi)/dx^2. Discretizing it with e.g. first order upwind and centered differencing yields
  //
  //    phi^(k+1)_i = phi^k_i - dt*v*[phi^k_i - phi^k_(i-1)]/dx + dt*D*[phi^k_(i+1) - 2*phi^k_i + phi^k_(i-1)]/(dx*dx).
  //
  // Setting phi_(i-1) = phi_(i+1) = 0 yields phi^(k+1)_i = phi^k_i - (dt*v/dx)*phi^k_i - 2*D*dt/(dx*dx)*phi_i^k. This is positive definite only if
  //
  //    (1 - dt*v/dx - 2*D*dt/(dx*dx)) > 0.0
  //
  // which yields
  //
  //    dt <= 1/(v/dx + 2*D/(dx*dx).
  //
  // Generalization to 3D yields dt <= 1/[(|vx|+|vy|+|vz)/dx + 2*d*D/(dx*dx)]. The below kernels compute the factor in all grid cells, and we return the
  // smallest value.

  Real minDt = std::numeric_limits<Real>::max();

  // Default to advection or diffusion time steps if the solver is only advective/diffusive.
  if (m_isMobile && !m_isDiffusive) {
    minDt = this->computeAdvectionDt();
  }
  else if (!m_isMobile && m_isDiffusive) {
    minDt = this->computeDiffusionDt();
  }
  else if (m_isMobile && m_isDiffusive) {

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const Real               dx    = m_amr->getDx()[lvl];
      const Real               dx2   = dx * dx;
      const DataIterator&      dit   = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(min : minDt)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        const EBCellFAB& velo       = (*m_cellVelocity[lvl])[din];
        const EBFluxFAB& diffCoFace = (*m_faceCenteredDiffusionCoefficient[lvl])[din];
        const EBISBox&   ebisbox    = ebisl[din];
        const Box        cellBox    = dbl[din];

        VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[din];

        // Single-valued data.
        const BaseFab<Real>& veloReg = velo.getSingleValuedFAB();

        // Regular kernel.
        auto regularKernel = [&](const IntVect& iv) -> void {
          if (ebisbox.isRegular(iv)) {
            Real v = 0.0;
            Real D = 0.0;

            // Compute |v|
            for (int dir = 0; dir < SpaceDim; dir++) {
              v += std::abs(veloReg(iv, dir));
            }

            // Get biggest diffusion coefficient in this cell.
            for (int dir = 0; dir < SpaceDim; dir++) {
              D = std::max(D, diffCoFace[dir].getSingleValuedFAB()(iv, m_comp));
              D = std::max(D, diffCoFace[dir].getSingleValuedFAB()(iv + BASISV(dir), m_comp));
            }

            const Real curDt = 1. / (v / dx + 2 * D * SpaceDim / dx2);

            minDt = std::min(minDt, curDt);
          }
        };

        // Irregular kernel.
        auto irregularKernel = [&](const VolIndex& vof) -> void {
          Real v = 0.0;
          Real D = 0.0;

          // Compute |v|
          for (int dir = 0; dir < SpaceDim; dir++) {
            v += std::abs(velo(vof, dir));
          }

          // Compute largest D
          for (int dir = 0; dir < SpaceDim; dir++) {
            const EBFaceFAB& Dface = diffCoFace[dir];

            for (SideIterator sit; sit.ok(); ++sit) {
              const std::vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, sit()).stdVector();

              for (const auto& face : faces) {
                D = std::max(D, Dface(face, m_comp));
              }
            }
          }

          const Real curDt = 1. / (v / dx + 2 * D * SpaceDim / dx2);

          minDt = std::min(minDt, curDt);
        };

        // Execute the kernels.
        BoxLoops::loop(dbl[din], regularKernel);
        BoxLoops::loop(vofit, irregularKernel);
      }
    }
  }

  return ParallelOps::min(minDt);
}

Real
CdrSolver::computeSourceDt(const Real a_max, const Real a_tolerance)
{
  CH_TIME("CdrSolver::computeSourceDt(Real, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeSourceDt(Real, Real)" << endl;
  }

  // This routine computes the time step dt = phi/source but only for cells where phi lies within a_tolerance*a_max. It's an old
  // routine -- and one that is hardly ever used.

  Real minDt = std::numeric_limits<Real>::max();

  if (a_max > 0.0) {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(min : minDt)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        const EBCellFAB& phi     = (*m_phi[lvl])[din];
        const EBCellFAB& source  = (*m_source[lvl])[din];
        const Box        cellBox = dbl[din];

        VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[din];

        const BaseFab<Real>& phiReg = phi.getSingleValuedFAB();
        const BaseFab<Real>& srcReg = source.getSingleValuedFAB();

        // Regular kernel.
        auto regularKernel = [&](const IntVect& iv) -> void {
          const Real curPhi = phiReg(iv, m_comp);
          const Real curSrc = srcReg(iv, m_comp);

          if (curPhi > a_tolerance * a_max && curSrc > 0.0) {
            minDt = std::min(minDt, std::abs(curPhi / curSrc));
          }
        };

        // Irregular kernel.
        auto irregularKernel = [&](const VolIndex& vof) -> void {
          const Real curPhi = phi(vof, m_comp);
          const Real curSrc = source(vof, m_comp);

          if (curPhi > a_tolerance * a_max && curSrc > 0.0) {
            minDt = std::min(minDt, std::abs(curPhi / curSrc));
          }
        };

        // Execute kernels.
        BoxLoops::loop(dbl[din], regularKernel);
        BoxLoops::loop(vofit, irregularKernel);
      }
    }
  }

  return ParallelOps::min(minDt);
}

void
CdrSolver::weightedUpwind(EBAMRCellData& a_weightedUpwindPhi, const int a_pow)
{
  CH_TIME("CdrSolver::weightedUpwind()");
  if (m_verbosity > 5) {
    pout() << m_name + "::weightedUpwind()" << endl;
  }

  if (m_isMobile) {
    EBAMRCellData scratch;
    m_amr->allocate(scratch, m_realm, m_phase, 1);

    m_amr->conservativeAverage(m_cellVelocity, m_realm, m_phase);
    m_amr->interpGhost(m_cellVelocity, m_realm, m_phase);

    m_amr->conservativeAverage(m_phi, m_realm, m_phase);
    m_amr->interpGhost(m_phi, m_realm, m_phase);

    // Compute velocity on faces and EBs. The data is currently face-centered.
    this->averageVelocityToFaces(m_faceVelocity, m_cellVelocity);
    this->advectToFaces(m_faceStates, m_phi, 0.0);

    DataOps::setValue(scratch, 0.0); // Used to store sum(alpha*v)

    // Grid loop.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {

      const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const DataIterator&      dit   = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        const Box&     cellBox = dbl[din];
        const EBISBox& ebisBox = ebisl[din];
        VoFIterator&   vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[din];

        EBCellFAB&       sumPhi    = (*a_weightedUpwindPhi[lvl])[din];
        EBCellFAB&       sumWeight = (*scratch[lvl])[din];
        const EBFluxFAB& facePhi   = (*m_faceStates[lvl])[din];
        const EBFluxFAB& faceVel   = (*m_faceVelocity[lvl])[din];

        sumPhi.setVal(0.0);
        sumWeight.setVal(0.0);

        // Regular cells.
        BaseFab<Real>& regSumPhi    = sumPhi.getSingleValuedFAB();
        BaseFab<Real>& regSumWeight = sumWeight.getSingleValuedFAB();

        for (int dir = 0; dir < SpaceDim; dir++) {
          const BaseFab<Real>& regFacePhi = facePhi[dir].getSingleValuedFAB();
          const BaseFab<Real>& regFaceVel = faceVel[dir].getSingleValuedFAB();

          // Regular kernel. Note that we execute the kernel over the cell-centered box.
          auto regularKernel = [&](const IntVect& iv) -> void {
            const Real& phiLo = regFacePhi(iv, m_comp);
            const Real& phiHi = regFacePhi(iv + BASISV(dir), m_comp);

            const Real& velLo = regFaceVel(iv, m_comp);
            const Real& velHi = regFaceVel(iv + BASISV(dir), m_comp);

            if (velLo > 0.0) {
              regSumWeight(iv, m_comp) += std::abs(std::pow(velLo, a_pow));
              regSumPhi(iv, m_comp) += std::abs(std::pow(velLo, a_pow)) * phiLo;
            }

            if (velHi < 0.0) {
              regSumWeight(iv, m_comp) += std::abs(std::pow(velHi, a_pow));
              regSumPhi(iv, m_comp) += std::abs(std::pow(velHi, a_pow)) * phiHi;
            }
          };

          // Execute the kernel
          BoxLoops::loop(cellBox, regularKernel);
        }

        // Irregular cells. This is a bit more involved. Note that we do want to compute everything at face centroids since
        // that is the flux that makes its way into the cells anyways when we advect. But, m_faceStates are states on the face
        // centers. So there's that.
        auto irregularKernel = [&](const VolIndex& vof) -> void {
          sumPhi(vof, m_comp)    = 0.0;
          sumWeight(vof, m_comp) = 0.0;

          for (int dir = 0; dir < SpaceDim; dir++) {
            const Vector<FaceIndex> facesLo = ebisBox.getFaces(vof, dir, Side::Lo);
            const Vector<FaceIndex> facesHi = ebisBox.getFaces(vof, dir, Side::Hi);

            // Add contribution from cut-cell faces in the low side. Observe that if the upwind value of phi is outside the domain
            // then the "upwinded" value is just the cell-centered value.
            for (int iface = 0; iface < facesLo.size(); iface++) {
              const FaceIndex& faceLo   = facesLo[iface];
              const VolIndex&  vofLo    = faceLo.getVoF(Side::Lo);
              const Real       areaFrac = ebisBox.areaFrac(faceLo);

              const FaceStencil& interpSten = (*m_interpStencils[dir][lvl])[din](faceLo, m_comp);

              Real phiLo = 0.0;
              Real velLo = 0.0;

              for (int i = 0; i < interpSten.size(); i++) {
                const FaceIndex& interpFace   = interpSten.face(i);
                const Real&      interpWeight = interpSten.weight(i);

                phiLo += interpWeight * facePhi[dir](interpFace, m_comp);
                velLo += interpWeight * faceVel[dir](interpFace, m_comp);
              }

              if (velLo > 0.0) {
                sumWeight(vof, m_comp) += areaFrac * std::abs(std::pow(velLo, a_pow));
                sumPhi(vof, m_comp) += areaFrac * std::abs(std::pow(velLo, a_pow)) * phiLo;
              }
            }

            // Add contribution from cut-cell faces on the high side.
            for (int iface = 0; iface < facesHi.size(); iface++) {
              const FaceIndex& faceHi   = facesHi[iface];
              const VolIndex&  vofHi    = faceHi.getVoF(Side::Hi);
              const Real       areaFrac = ebisBox.areaFrac(faceHi);

              const FaceStencil& interpSten = (*m_interpStencils[dir][lvl])[din](faceHi, m_comp);

              Real phiHi = 0.0;
              Real velHi = 0.0;

              for (int i = 0; i < interpSten.size(); i++) {
                const FaceIndex& interpFace   = interpSten.face(i);
                const Real&      interpWeight = interpSten.weight(i);

                phiHi += interpWeight * facePhi[dir](interpFace, m_comp);
                velHi += interpWeight * faceVel[dir](interpFace, m_comp);
              }

              if (velHi < 0.0) {
                sumWeight(vof, m_comp) += areaFrac * std::abs(std::pow(velHi, a_pow));
                sumPhi(vof, m_comp) += areaFrac * std::abs(std::pow(velHi, a_pow)) * phiHi;
              }
            }
          }

          // If we don't have an inflow face set the upwind stuff to zero.
          if (!(sumWeight(vof, m_comp) > 0.0)) {
            sumPhi(vof, m_comp)    = 0.0;
            sumWeight(vof, m_comp) = 1.0;
          }
        };

        // Execute kernel.
        BoxLoops::loop(vofit, irregularKernel);
      }
    }

    // Divide. Set to zero m_phi if there are no inflow faces.
    EBAMRCellData zero;
    m_amr->allocate(zero, m_realm, m_phase, m_nComp);
    DataOps::setValue(zero, 0.0);
    DataOps::divideFallback(a_weightedUpwindPhi, scratch, zero);

    m_amr->conservativeAverage(a_weightedUpwindPhi, m_realm, m_phase);
    m_amr->interpGhost(a_weightedUpwindPhi, m_realm, m_phase);
  }
  else {
    DataOps::copy(a_weightedUpwindPhi, m_phi);
  }
}

Real
CdrSolver::computeMass()
{
  CH_TIME("CdrSolver::computeMass()");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeMass()" << endl;
  }

  return this->computeMass(m_phi);
}

Real
CdrSolver::computeMass(const EBAMRCellData& a_phi, const bool a_kappaScale)
{
  CH_TIME("CdrSolver::computeMass(EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeMass(EBAMRCellData)" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);

  Real mass = 0.0;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real               dx    = m_amr->getDx()[lvl];
    const Real               dxVol = std::pow(dx, SpaceDim);

    const DataIterator& dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(+ : mass)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box            cellbox    = dbl[din];
      const EBISBox&       ebisbox    = ebisl[din];
      const BaseFab<bool>& validCells = (*m_amr->getValidCells(m_realm)[lvl])[din];

      const EBCellFAB& phi    = (*a_phi[lvl])[din];
      const FArrayBox& phiReg = phi.getFArrayBox();

      // Kernel definitions.
      auto regularKernel = [&](const IntVect& iv) -> void {
        if (validCells(iv) && ebisbox.isRegular(iv)) {
          mass += phiReg(iv, 0) * dxVol;
        }
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const IntVect iv = vof.gridIndex();

        if (validCells(iv)) {
          const Real volFrac = a_kappaScale ? ebisbox.volFrac(vof) : 1.0;

          mass += phi(vof, 0) * dxVol * volFrac;
        }
      };

      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[din];

      BoxLoops::loop(cellbox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }

  return ParallelOps::sum(mass);
}

Real
CdrSolver::computeCharge()
{
  CH_TIME("CdrSolver::computeCharge()");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeCharge()" << endl;
  }

  // TLDR: Compute the mass on the coarsest level and multiply by charge number.

  const Real Q = this->computeMass() * m_species->getChargeNumber();

  return Q;
}

bool
CdrSolver::isDiffusive()
{
  CH_TIME("CdrSolver::isDiffusive()");
  if (m_verbosity > 5) {
    pout() << m_name + "::isDiffusive()" << endl;
  }

  return m_isDiffusive;
}

bool
CdrSolver::isMobile()
{
  CH_TIME("CdrSolver::isMobile()");
  if (m_verbosity > 5) {
    pout() << m_name + "::isMobile()" << endl;
  }

  return m_isMobile;
}

EBAMRCellData&
CdrSolver::getPhi()
{
  CH_TIME("CdrSolver::getPhi()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getPhi()" << endl;
  }

  return m_phi;
}

EBAMRCellData&
CdrSolver::getSource()
{
  CH_TIME("CdrSolver::getSource()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getSource()" << endl;
  }

  return m_source;
}

EBAMRCellData&
CdrSolver::getCellCenteredVelocity()
{
  CH_TIME("CdrSolver::getCellCenteredVelocity()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getCellCenteredVelocity()" << endl;
  }

  return m_cellVelocity;
}

EBAMRFluxData&
CdrSolver::getFaceCenteredVelocity()
{
  CH_TIME("CdrSolver::getFaceCenteredVelocity()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getFaceCenteredVelocity()" << endl;
  }

  return m_faceVelocity;
}

EBAMRIVData&
CdrSolver::getEbCenteredVelocity()
{
  CH_TIME("CdrSolver::getEbCenteredVelocity()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getEbCenteredVelocity()" << endl;
  }

  return m_ebVelocity;
}

EBAMRCellData&
CdrSolver::getCellCenteredDiffusionCoefficient()
{
  CH_TIME("CdrSolver::getCellCenteredDiffusionCoefficient()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getCellCenteredDiffusionCoefficient()" << endl;
  }

  return m_cellCenteredDiffusionCoefficient;
}

EBAMRFluxData&
CdrSolver::getFaceCenteredDiffusionCoefficient()
{
  CH_TIME("CdrSolver::getFaceCenteredDiffusionCoefficient()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getFaceCenteredDiffusionCoefficient()" << endl;
  }

  return m_faceCenteredDiffusionCoefficient;
}

EBAMRIVData&
CdrSolver::getEbCenteredDiffusionCoefficient()
{
  CH_TIME("CdrSolver::getEbCenteredDiffusionCoefficient()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getEbCenteredDiffusionCoefficient()" << endl;
  }

  return m_ebCenteredDiffusionCoefficient;
}

EBAMRIVData&
CdrSolver::getEbFlux()
{
  CH_TIME("CdrSolver::getEbFlux()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getEbFlux()" << endl;
  }

  return m_ebFlux;
}

EBAMRIFData&
CdrSolver::getDomainFlux()
{
  CH_TIME("CdrSolver::getDomainFlux()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getDomainFlux()" << endl;
  }

  return m_domainFlux;
}

void
CdrSolver::extrapolateAdvectiveFluxToEB() noexcept
{
  CH_TIME("CdrSolver::extrapolateAdvectiveFluxToEB()");
  if (m_verbosity > 5) {
    pout() << m_name + "::extrapolateAdvectiveFluxToEB()" << endl;
  }

  this->extrapolateAdvectiveFluxToEB(m_ebFlux);
}

void
CdrSolver::extrapolateAdvectiveFluxToEB(EBAMRIVData& a_ebFlux) const noexcept
{
  CH_TIME("CdrSolver::extrapolateAdvectiveFluxToEB(EBAMRIVData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::extrapolateadvectivefluxtoeb(EBAMRIVData)" << endl;
  }

  DataOps::setValue(a_ebFlux, 0.0);

  if (m_isMobile) {

    // Compute v * phi on cell centers.
    EBAMRCellData flux;
    EBAMRIVData   fluxEB;

    m_amr->allocate(flux, m_realm, m_phase, SpaceDim);
    m_amr->allocate(fluxEB, m_realm, m_phase, SpaceDim);

    DataOps::copy(flux, m_cellVelocity);
    DataOps::multiplyScalar(flux, m_phi);

    m_amr->conservativeAverage(flux, m_realm, m_phase);
    m_amr->interpGhost(flux, m_realm, m_phase);

    // Extrapolate to the EB.
    m_amr->interpToEB(fluxEB, flux, m_realm, m_phase);

    // Project along normal vector.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      CH_assert(a_ebFlux[lvl]->nComp() == 1);

      const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const DataIterator&      dit   = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        const Box      cellBox = dbl[din];
        const EBISBox& ebisBox = ebisl[din];

        BaseIVFAB<Real>&       scalarFlux = (*a_ebFlux[lvl])[din];
        const BaseIVFAB<Real>& vectorFlux = (*fluxEB[lvl])[din];

        auto irregularKernel = [&](const VolIndex& vof) {
          const RealVect& normal = ebisBox.normal(vof);

          for (int dir = 0; dir < SpaceDim; dir++) {
            scalarFlux(vof, 0) -= vectorFlux(vof, dir) * normal[dir];
          }
        };

        VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[din];

        BoxLoops::loop(vofit, irregularKernel);
      }
    }
  }
}

std::string
CdrSolver::makeBcString(const int a_dir, const Side::LoHiSide a_side) const
{
  CH_TIME("CdrSolver::makeBcString(int, Side::LoHiSide)");
  if (m_verbosity > 5) {
    pout() << "CdrSolver::makeBcString(int, Side::LoHiSide)" << endl;
  }

  std::string strDir;
  std::string strSide;

  if (a_dir == 0) {
    strDir = "x";
  }
  else if (a_dir == 1) {
    strDir = "y";
  }
  else if (a_dir == 2) {
    strDir = "z";
  }

  if (a_side == Side::Lo) {
    strSide = "lo";
  }
  else if (a_side == Side::Hi) {
    strSide = "hi";
  }

  const std::string ret = std::string("bc.") + strDir + std::string(".") + strSide;

  return ret;
}

void
CdrSolver::parseDomainBc()
{
  CH_TIME("CdrSolver::parseDomainBc()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseDomainBc()" << endl;
  }

  ParmParse pp(m_className.c_str());

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      const Side::LoHiSide side = sit();

      const CdrDomainBC::DomainSide domainSide = std::make_pair(dir, side);
      const std::string             bcString   = this->makeBcString(dir, side);

      std::string str;
      pp.get(bcString.c_str(), str);

      // Set domain bc from input script evaluations.
      if (str == "data") {
        this->setDomainBcType(domainSide, CdrDomainBC::BcType::DataBased);
      }
      else if (str == "function") {
        this->setDomainBcType(domainSide, CdrDomainBC::BcType::Function);
      }
      else if (str == "wall") {
        this->setDomainBcType(domainSide, CdrDomainBC::BcType::Wall);
      }
      else if (str == "outflow") {
        this->setDomainBcType(domainSide, CdrDomainBC::BcType::Outflow);
      }
      else if (str == "solver") {
        this->setDomainBcType(domainSide, CdrDomainBC::BcType::Solver);
      }
      else {
        const std::string errorString = "CdrSolver::parseDomain BC - bad argument '" + bcString + "'";
        MayDay::Error(errorString.c_str());
      }
    }
  }
}

void
CdrSolver::parseDivergenceComputation()
{
  CH_TIME("CdrSolver::parseDivergenceComputation()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseDivergenceComputation()" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;

  // Get redistribution type.
  pp.get("which_redistribution", str);
  if (str == "volume") {
    m_whichRedistribution = Redistribution::VolumeWeighted;
  }
  else if (str == "none") {
    m_whichRedistribution = Redistribution::None;
  }
  else {
    MayDay::Error("CdrSolver::parseDivergenceComputation -- logic bust. Must specify 'volume' or 'none'");
  }

  // Blend conservation or not.
  pp.get("blend_conservation", m_blendConservation);
}

void
CdrSolver::parsePlotVariables()
{
  CH_TIME("CdrSolver::parsePlotVariables()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parsePlotVariables()" << endl;
  }

  m_plotPhi                  = false;
  m_plotVelocity             = false;
  m_plotDiffusionCoefficient = false;
  m_plotSource               = false;
  m_plotEbFlux               = false;

  // Parse plot variables from command line
  ParmParse pp(m_className.c_str());
  const int num = pp.countval("plt_vars");

  if (num > 0) {
    Vector<std::string> str(num);
    pp.getarr("plt_vars", str, 0, num);

    // Set plot variables
    for (int i = 0; i < num; i++) {
      if (str[i] == "phi") {
        m_plotPhi = true;
      }
      else if (str[i] == "vel") {
        m_plotVelocity = true;
      }
      else if (str[i] == "dco") {
        m_plotDiffusionCoefficient = true;
      }
      else if (str[i] == "src") {
        m_plotSource = true;
      }
      else if (str[i] == "ebflux") {
        m_plotEbFlux = true;
      }
    }
  }
}

void
CdrSolver::gwnDiffusionSource(EBAMRCellData& a_noiseSource, const EBAMRCellData& a_cellPhi)
{
  CH_TIME("CdrSolver::gwnDiffusionSource(EBAMRCellData, EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::gwnDiffusionSource(EBAMRCellData, EBAMRCellData)" << endl;
  }

  CH_assert(a_noiseSource[0]->nComp() == 1);
  CH_assert(a_cellPhi[0]->nComp() == 1);

  // This routine computes a Gaussian white noise source term for the fluctuating convection-diffusion equation. The source
  // term enters the convection-diffusion equation as
  //
  //    d(phi)/dt = Div(D*Grad(phi) - Z*sqrt(2*D*phi)) where Z is a Gaussian random field.
  //
  // In this routine we compute the finite volume approximation to Div(Z*sqrt(2*D*phi)) using standard finite volume routines. The only difference
  // is that the flux is defined as Z*sqrt(2*D*phi) but we don't have phi on face centers. Worse, setting phi to be the arithmetic average can easily
  // lead to negative density, so we use the Heaviside smoothing from Kim et. al. The paper is
  // "Stochastic Simulation of Reaction-Diffusion Systems: A Fluctuating-Hydrodynamics Approach", 10.1063/1.4978775.
  //
  // This does not guarantee against negative densities, but its surely better than the naive approach. Note that once we have the fluctuating term, we simply
  // compute Div(G).

  if (m_isDiffusive) {
    EBAMRFluxData scratchFluxOne;
    EBAMRFluxData scratchFluxTwo;

    m_amr->allocate(scratchFluxOne, m_realm, m_phase, m_nComp);
    m_amr->allocate(scratchFluxTwo, m_realm, m_phase, m_nComp);

    // scratchFluxOne = phis on faces (smoothing as to avoid negative densities)
    this->smoothHeavisideFaces(scratchFluxOne, a_cellPhi);
    DataOps::multiply(scratchFluxOne, m_faceCenteredDiffusionCoefficient); // scratchFluxOne = D*phis
    DataOps::scale(scratchFluxOne, 2.0);                                   // scratchFluxOne = 2*D*phis
    DataOps::squareRoot(scratchFluxOne);                                   // scratchFluxOne = sqrt(2*D*phis)

#ifndef NDEBUG
    Real max;
    Real min;

    DataOps::getMaxMin(max, min, scratchFluxOne, 0);

    if (min < 0.0 || max < 0.0) {
      MayDay::Abort("CdrSolver::gwnDiffusionSource - negative face value");
    }
#endif

    // Let scratchFluxOne holds the fluctuating cell-centered flux Z*sqrt(2*D*phi) and then compute the finite volume approximation of the
    // divergence
    this->fillGwn(scratchFluxTwo, 1.0);
    DataOps::multiply(scratchFluxOne, scratchFluxTwo);
    this->computeDivG(a_noiseSource, scratchFluxOne, m_ebZero, false);
  }
  else {
    DataOps::setValue(a_noiseSource, 0.0);
  }
}

void
CdrSolver::smoothHeavisideFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_cellPhi)
{
  CH_TIME("CdrSolver::smoothHeavisideFaces(EBAMRFluxData, EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::smoothHeavisideFaces(EBAMRFluxData, EBAMRCellData)" << endl;
  }

  CH_assert(a_facePhi[0]->nComp() == 1);
  CH_assert(a_cellPhi[0]->nComp() == 1);

  // This routine is taken from  Kim et. al. "Stochastic Simulation of Reaction-Diffusion Systems: A Fluctuating-Hydrodynamics Approach" (10.1063/1.4978775).
  // It is inspired by a desire to gradually turn off fluctuations as the number of particles in a grid become small, as to avoid negative densities. So, we
  // compute the value of phi on faces with the following rules:
  //
  //    1. If theres more than one particle in the cells, we take the arithmetic average as usual.
  //    2. If there's between zero and one particle in the cells, we use an averaging function phis = 0.5*(phiLo + phiHi) * loFactor * hiFactor
  //       where loFactor and hiFactor are the number of particles in the grid cell.

  // Loop over levels.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real               dx    = m_amr->getDx()[lvl];
    const Real               vol   = pow(dx, SpaceDim);
    const DataIterator&      dit   = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box&        cellBox  = dbl[din];
      const EBISBox&    ebisbox  = ebisl[din];
      const EBGraph&    ebgraph  = ebisbox.getEBGraph();
      const IntVectSet& irregIVS = ebisbox.getIrregIVS(cellBox);

      for (int dir = 0; dir < SpaceDim; dir++) {
        EBFaceFAB&       facePhi = (*a_facePhi[lvl])[din][dir];
        const EBCellFAB& cellPhi = (*a_cellPhi[lvl])[din];

        BaseFab<Real>&       faceReg = facePhi.getSingleValuedFAB();
        const BaseFab<Real>& cellReg = cellPhi.getSingleValuedFAB();

        // Regular kernel
        auto regularKernel = [&](const IntVect& iv) -> void {
          const Real& phiLo = cellReg(iv - BASISV(dir), m_comp);
          const Real& phiHi = cellReg(iv, m_comp);

          const Real loVal = vol * phiLo;
          const Real hiVal = vol * phiHi;

          Real Hlo = 0.0;
          Real Hhi = 0.0;

          if (loVal <= 0.0) {
            Hlo = 0.0;
          }
          else if (loVal >= 1.0) {
            Hlo = 1.0;
          }
          else {
            Hlo = loVal;
          }

          if (hiVal <= 0.0) {
            Hhi = 0.0;
          }
          else if (hiVal >= 1.0) {
            Hhi = 1.0;
          }
          else {
            Hhi = hiVal;
          }

          faceReg(iv, m_comp) = 0.5 * (phiLo + phiHi) * Hlo * Hhi;
        };

        // Irregular kernel. Same as the above but we need to fetch faces explicitly.
        auto irregularKernel = [&](const FaceIndex& face) -> void {
          const VolIndex loVof = face.getVoF(Side::Lo);
          const VolIndex hiVof = face.getVoF(Side::Hi);

          const Real loVal = std::max((Real)0.0, cellPhi(loVof, m_comp));
          const Real hiVal = std::max((Real)0.0, cellPhi(hiVof, m_comp));

          Real Hlo = 0.0;
          Real Hhi = 0.0;

          if (loVal * vol <= 0.0) {
            Hlo = 0.0;
          }
          else if (loVal * vol >= 1.0) {
            Hlo = 1.0;
          }
          else {
            Hlo = loVal * vol;
          }

          if (hiVal * vol <= 0.0) {
            Hhi = 0.0;
          }
          else if (hiVal * vol >= 1.0) {
            Hhi = 1.0;
          }
          else {
            Hhi = hiVal * vol;
          }

          facePhi(face, m_comp) = 0.5 * (hiVal + loVal) * Hlo * Hhi;
        };

        // These are the computation regions for the kernels.
        const Box    faceBox = surroundingNodes(cellBox, dir);
        FaceIterator faceit(irregIVS, ebgraph, dir, FaceStop::SurroundingNoBoundary);

        // Execute the kernels.
        BoxLoops::loop(faceBox, regularKernel);
        BoxLoops::loop(faceit, irregularKernel);
      }
    }

    // Covered faces are bogus.
    DataOps::setCoveredValue(*a_facePhi[lvl], 0.0);
  }

  // No random flux on domain faces.
  this->resetDomainFlux(a_facePhi);
}

void
CdrSolver::fillGwn(EBAMRFluxData& a_noise, const Real a_sigma)
{
  CH_TIME("CdrSolver::fillGwn(EBAMRFluxData, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::fillGwn(EBAMRFluxData, Real)" << endl;
  }

  CH_assert(a_noise[0]->nComp() == 1);

  // TLDR: This routine draws a random number from a Gaussian distribution at each cell face. This is used
  //       in FHD routines where we want to compute a finite volume approximation to the FHD diffusion noise
  //       term.

  // Gaussian white noise distribution -- this should be centered at 0, and we will usually but not necessarily have a_sigma=1
  std::normal_distribution<double> whiteNoise(0.0, a_sigma);

  // Initialize.
  DataOps::setValue(a_noise, 0.0);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real               dx    = m_amr->getDx()[lvl];
    const Real               vol   = pow(dx, SpaceDim);
    const Real               ivol  = sqrt(1. / vol);
    const DataIterator&      dit   = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box&        cellBox = dbl[din];
      const EBISBox&    ebisbox = ebisl[din];
      const EBGraph&    ebgraph = ebisbox.getEBGraph();
      const IntVectSet& irreg   = ebisbox.getIrregIVS(cellBox);

      for (int dir = 0; dir < SpaceDim; dir++) {
        EBFaceFAB&     noise    = (*a_noise[lvl])[din][dir];
        BaseFab<Real>& noiseReg = noise.getSingleValuedFAB();

        // Regular faces
        const Box    facebox = surroundingNodes(cellBox, dir);
        FaceIterator faceit(irreg, ebgraph, dir, FaceStop::SurroundingNoBoundary);

        // Regular kernel
        auto regularKernel = [&](const IntVect& iv) -> void {
          noiseReg(iv, m_comp) = Random::get(whiteNoise) * ivol;
        };

        // Irregular kernel
        auto irregularKernel = [&](const FaceIndex& face) -> void {
          noise(face, m_comp) = Random::get(whiteNoise) * ivol;
        };

        // Execute the kernels.
        BoxLoops::loop(facebox, regularKernel);
        BoxLoops::loop(faceit, irregularKernel);
      }
    }
  }
}

void
CdrSolver::parsePlotMode()
{
  CH_TIME("CdrSolver::parsePlotMode()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parsePlotMode()" << endl;
  }

  // Parses plot mode.

  ParmParse pp(m_className.c_str());

  m_plotNumbers = false;

  std::string str;
  pp.get("plot_mode", str);
  if (str == "density") {
    m_plotNumbers = false;
  }
  else if (str == "numbers") {
    m_plotNumbers = true;
  }
}

void
CdrSolver::parseRegridSlopes()
{
  CH_TIME("CdrSolver::parseRegridSlopes()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRegridSlopes()" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("use_regrid_slopes", m_regridSlopes);
}

#include <CD_NamespaceFooter.H>
