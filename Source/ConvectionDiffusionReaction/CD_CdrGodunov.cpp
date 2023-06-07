/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrGodunov.cpp
  @brief  Implementation of CD_CdrGodunov.H
  @author Robert Marskar
*/

// Chombo includes
#include <ExtrapAdvectBC.H>
#include <EBArith.H>
#include <ParmParse.H>

// Our includes
#include <CD_CdrGodunov.H>
#include <CD_DataOps.H>
#include <CD_ParallelOps.H>
#include "CD_NamespaceHeader.H"

CdrGodunov::CdrGodunov() : CdrMultigrid()
{
  CH_TIME("CdrGodunov::CdrGodunov()");

  // Class name and instantiation name.
  m_className = "CdrGodunov";
  m_name      = "CdrGodunov";
}

CdrGodunov::~CdrGodunov() { CH_TIME("CdrGodunov::~CdrGodunov()"); }

void
CdrGodunov::parseOptions()
{
  CH_TIME("CdrGodunov::parseOptions()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseOptions()" << endl;
  }

  this->parseDomainBc();              // Parses domain BC options
  this->parseSlopeLimiter();          // Parses slope limiter settings
  this->parsePlotVariables();         // Parses plot variables
  this->parsePlotMode();              // Parses plot mode
  this->parseMultigridSettings();     // Parses multigrid settings
  this->parseExtrapolateSourceTerm(); // Parses source term extrapolation for Godunov time extrapolation.
  this->parseDivergenceComputation(); // Parses non-conservative divergence blending
  this->parseRegridSlopes();          // Parses regrid slopes
}

void
CdrGodunov::parseRuntimeOptions()
{
  CH_TIME("CdrGodunov::parseRuntimeOptions()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRuntimeOptions()" << endl;
  }

  this->parseSlopeLimiter();          // Parses slope limiter (on/off)
  this->parsePlotVariables();         // Parses plot variables
  this->parsePlotMode();              // Parses plot mode
  this->parseMultigridSettings();     // Parses multigrid settings
  this->parseDomainBc();              // Parses domain BCs
  this->parseExtrapolateSourceTerm(); // Parses source term extrapolation
  this->parseDivergenceComputation(); // Parses non-conservative divergence blending.
  this->parseRegridSlopes();          // Parses regrid slopes
}

Real
CdrGodunov::computeAdvectionDt()
{
  CH_TIME("CdrGodunov::computeAdvectionDt()");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeAdvectionDt()" << endl;
  }

  // TLDR: For advection, Bell, Collela, and Glaz says we must have dt <= dx/max(|vx|, |vy|, |vz|). See these two papers for details:
  //
  //       Bell, Colella, Glaz, J. Comp. Phys 85 (257), 1989
  //       Minion, J. Comp. Phys 123 (435), 1996

  Real minDt = std::numeric_limits<Real>::max();

  if (m_isMobile) {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const Real               dx    = m_amr->getDx()[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        const Box        cellBox = dbl[dit()];
        const EBCellFAB& velo    = (*m_cellVelocity[lvl])[dit()];
        const EBISBox&   ebisBox = ebisl[dit()];

        VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

        // Regular grid data.
        const BaseFab<Real>& veloReg = velo.getSingleValuedFAB();

        // Compute dt = dx/(|vx|+|vy|+|vz|) and check if it's smaller than the smallest so far.
        auto regularKernel = [&](const IntVect& iv) -> void {
          Real velMax = 0.0;
          if (ebisBox.isRegular(iv)) {
            for (int dir = 0; dir < SpaceDim; dir++) {
              velMax = std::max(velMax, std::abs(veloReg(iv, dir)));
            }
          }

          if (velMax > 0.0) {
            minDt = std::min(dx / velMax, minDt);
          }
        };

        // Same kernel, but for cut-cells.
        auto irregularKernel = [&](const VolIndex& vof) -> void {
          Real velMax = 0.0;
          for (int dir = 0; dir < SpaceDim; dir++) {
            velMax = std::max(velMax, std::abs(velo(vof, dir)));
          }

          if (velMax > 0.0) {
            minDt = std::min(dx / velMax, minDt);
          }
        };

        // Execute the kernels.
        BoxLoops::loop(cellBox, regularKernel);
        BoxLoops::loop(vofit, irregularKernel);
      }
    }

    // If we are using MPI then ranks need to know of each other's time steps.
    minDt = ParallelOps::min(minDt);
  }

  return minDt;
}

void
CdrGodunov::parseExtrapolateSourceTerm()
{
  CH_TIME("CdrGodunov::parseExtrapolateSourceTerm()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseExtrapolateSourceTerm()" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("extrap_source", m_extrapolateSourceTerm);
}

void
CdrGodunov::parseSlopeLimiter()
{
  CH_TIME("CdrGodunov::parseSlopeLimiter()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseSlopeLimiter()" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("limit_slopes", m_limitSlopes);
}

void
CdrGodunov::allocate()
{
  CH_TIME("CdrSolver::allocate()");
  if (m_verbosity > 5) {
    pout() << m_name + "::allocate()" << endl;
  }

  // CdrMultigrid allocates everything except storage needed for the advection object.
  CdrMultigrid::allocate();

  // Allocate levelAdvect only if the solver is mobile. See Chombo design docs for how the EBAdvectLevelIntegrator operates.
  if (m_isMobile) {
    const Vector<RefCountedPtr<EBLevelGrid>>& eblgs       = m_amr->getEBLevelGrid(m_realm, m_phase);
    const Vector<int>&                        refRatios   = m_amr->getRefinementRatios();
    const Vector<Real>&                       dx          = m_amr->getDx();
    const int                                 finestLevel = m_amr->getFinestLevel();

    m_levelAdvect.resize(1 + finestLevel);

    for (int lvl = 0; lvl <= finestLevel; lvl++) {

      const bool hasCoar = lvl > 0;
      const bool hasFine = lvl < finestLevel;

      int         refRat = 1;
      EBLevelGrid coarEblg;
      if (hasCoar) {
        coarEblg = *eblgs[lvl - 1];
        refRat   = refRatios[lvl - 1];
      }

      // Note: There is a "bug" in the function signature in Chombo. The second-to-last argument is the slope limiter and not the EBCF things.
      const EBIndexSpace* const ebis = eblgs[lvl]->getEBIS();

      m_levelAdvect[lvl] = RefCountedPtr<EBAdvectLevelIntegrator>(new EBAdvectLevelIntegrator(*eblgs[lvl],
                                                                                              coarEblg,
                                                                                              refRat,
                                                                                              dx[lvl] * RealVect::Unit,
                                                                                              hasCoar,
                                                                                              hasFine,
                                                                                              false,
                                                                                              m_limitSlopes,
                                                                                              ebis));
    }
  }
}

void
CdrGodunov::advectToFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_cellPhi, const Real a_extrapDt)
{
  CH_TIME("CdrGodunov::advectToFaces(EBAMRFluxDat, EBAMRCellData, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::advectToFaces(EBAMRFluxDat, EBAMRCellData, Real)" << endl;
  }

  CH_assert(a_facePhi[0]->nComp() == 1);
  CH_assert(a_cellPhi[0]->nComp() == 1);

  // If we are extrapolating in time the source term will yield different states on the face centers. The source term is the source
  // S^k + Div(D*Grad(Phi)), which we add to the solver below.
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, 1);

  DataOps::setValue(scratch, 0.0);

  // Compute viscous source term for the advection.
#if 0
  if(m_isDiffusive && a_extrapDt > 0.0) {
    this->setupDiffusionSolver();

    this->computeKappaLphi(scratch, a_cellPhi);
  }
#endif

  // If user asks for it, add the source term to the extrapolation.
  if (m_extrapolateSourceTerm && a_extrapDt > 0.0) {
    DataOps::incr(scratch, m_source, 1.0);
  }

  m_amr->conservativeAverage(scratch, m_realm, m_phase);
  m_amr->interpGhost(scratch, m_realm, m_phase);

  // This code extrapolates the cell-centered state to face centers on every grid level, in both space and time.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      EBFluxFAB&       facePhi = (*a_facePhi[lvl])[dit()];
      const EBCellFAB& cellPhi = (*a_cellPhi[lvl])[dit()];
      const EBCellFAB& cellVel = (*m_cellVelocity[lvl])[dit()];
      const EBFluxFAB& faceVel = (*m_faceVelocity[lvl])[dit()];
      const EBCellFAB& source  = (*scratch[lvl])[dit()];
      const Real       time    = 0.0;

      EBAdvectPatchIntegrator& ebAdvectPatch = m_levelAdvect[lvl]->getPatchAdvect(dit());

      // These are settings for EBAdvectPatchIntegrator -- it's not a very pretty design but the object has settings
      // that permits it to run advection code (through setDoingVel(0)).
      ebAdvectPatch.setVelocities(cellVel, faceVel);
      ebAdvectPatch.setDoingVel(0);
      ebAdvectPatch.setCurComp(m_comp);
      ebAdvectPatch.setEBPhysIBC(ExtrapAdvectBCFactory());

      // Extrapolate to face-centers. The face-centered states are Godunov-style extrapolated in time to a_extrapDt.
      ebAdvectPatch.extrapolateBCG(facePhi, cellPhi, source, dit(), time, a_extrapDt);
    }
  }
}

#include <CD_NamespaceFooter.H>
