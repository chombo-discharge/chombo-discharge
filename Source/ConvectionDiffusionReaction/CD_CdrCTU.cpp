/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_CdrCTU.cpp
 * @brief  Implementation of CD_CdrCTU.H
 * @author Robert Marskar
 */

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_CdrCTU.H>
#include <CD_BoxLoops.H>
#include <CD_DataOps.H>
#include <CD_ParallelOps.H>
#include <CD_NamespaceHeader.H>

CdrCTU::CdrCTU()
  : m_limiter(Limiter::MonotonizedCentral), m_useCTU(true), m_ebInflowState(false), m_fluxCorrection(false)
{
  CH_TIME("CdrCTU::CdrCTU()");

  // Class and object name
  m_className = "CdrCTU";
  m_name      = "CdrCTU";
}

CdrCTU::~CdrCTU()
{
  CH_TIME("CdrCTU::~CdrCTU()");
}

void
CdrCTU::parseOptions()
{
  CH_TIME("CdrCTU::parseOptions()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseOptions()" << endl;
  }

  this->parsePlotMode();              // Parses plot mode
  this->parseDomainBc();              // Parses domain BC options
  this->parseSlopeLimiter();          // Parses slope limiter settings
  this->parsePlotVariables();         // Parses plot variables
  this->parseMultigridSettings();     // Parses multigrid settings.
  this->parseDivergenceComputation(); // Non-conservative divergence blending
  this->parseRegridSlopes();          // Parses regrid slopes
}

void
CdrCTU::parseRuntimeOptions()
{
  CH_TIME("CdrCTU::parseRuntimeOptions()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRuntimeOptions()" << endl;
  }

  this->parsePlotMode();              // Parses plot mode
  this->parseDomainBc();              // Parses domain BC options
  this->parseSlopeLimiter();          // Parses slope limiter settings
  this->parsePlotVariables();         // Parses plot variables
  this->parseMultigridSettings();     // Parses multigrid settings.
  this->parseDivergenceComputation(); // Non-conservative divergence blending.
  this->parseRegridSlopes();          // Parses regrid slopes
}

void
CdrCTU::allocate()
{
  CH_TIME("CdrCTU::allocate()");
  if (m_verbosity > 5) {
    pout() << m_name + "::allocate()" << endl;
  }

  CdrMultigrid::allocate();
  this->defineIterators();
}

void
CdrCTU::preRegrid(const int a_lbase, const int a_oldFinestLevel)
{
  CH_TIME("CdrCTU::preRegrid(int, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::preRegrid(int, int)" << endl;
  }

  CdrMultigrid::preRegrid(a_lbase, a_oldFinestLevel);
  m_grownVofIter.resize(0);
  m_levelCells.resize(0);
}

void
CdrCTU::regrid(const int a_lbase, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("CdrCTU::regrid(int, int, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::regrid(int, int, int)" << endl;
  }

  CdrMultigrid::regrid(a_lbase, a_oldFinestLevel, a_newFinestLevel);
  this->defineIterators();
}

void
CdrCTU::defineIterators()
{
  CH_TIME("CdrCTU::defineIterators()");
  if (m_verbosity > 5) {
    pout() << m_name + "::defineIterators()" << endl;
  }

  const int finest = m_amr->getFinestLevel();

  m_grownVofIter.resize(1 + finest);
  m_levelCells.resize(1 + finest);

  for (int lvl = 0; lvl <= finest; lvl++) {
    const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const ProblemDomain&     domain = m_amr->getDomains()[lvl];

    m_grownVofIter[lvl] = RefCountedPtr<LayoutData<VoFIterator>>(new LayoutData<VoFIterator>(dbl));

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      const Box      grownBox = grow(dbl[dit], 1) & domain;
      const EBISBox& ebisbox  = ebisl[dit()];
      (*m_grownVofIter[lvl])[dit()].define(ebisbox.getIrregIVS(grownBox), ebisbox.getEBGraph());
    }

    // 1 in the cells of this level's grids, ghost cells covered by other boxes included, and 0 elsewhere.
    m_levelCells[lvl] = RefCountedPtr<LevelData<FArrayBox>>(new LevelData<FArrayBox>(dbl, 1, IntVect::Unit));

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      (*m_levelCells[lvl])[dit()].setVal(0.0);
      (*m_levelCells[lvl])[dit()].setVal(1.0, dbl[dit()], 0, 1);
    }

    m_levelCells[lvl]->exchange();
  }
}

Real
CdrCTU::computeAdvectionDt()
{
  CH_TIME("CdrCTU::computeAdvectionDt()");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeAdvectionDt()" << endl;
  }

  Real minDt = std::numeric_limits<Real>::max();

  if (!m_useCTU) {
    minDt = CdrMultigrid::computeAdvectionDt();
  }
  else {

    // TLDR: For advection, Bell, Collela, and Glaz says we must have dt <= dx/max(|vx|, |vy|, |vz|). See these three
    // papers for details:
    //
    //       Colella, J. Comp. Phys. 87 (171-200), 1990
    //       Bell, Colella, Glaz, J. Comp. Phys 85 (257), 1989
    //       Minion, J. Comp. Phys 123 (435), 1996
    if (m_isMobile) {
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
        const Real               dx  = m_amr->getDx()[lvl];
        const DataIterator&      dit = dbl.dataIterator();

        const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(min : minDt)
        for (int mybox = 0; mybox < nbox; mybox++) {
          const DataIndex& din = dit[mybox];

          const Box        cellBox = dbl[din];
          const EBCellFAB& velo    = (*m_cellVelocity[lvl])[din];

          // Use the precomputed not-covered mask (1 in regular/irregular cells, 0 in covered) instead of a per-cell
          // out-of-line EBISBox::isRegular query, and iterate only the MULTI-valued cut-cells below. Singly-cut cells
          // are handled by the regular box loop (their single-valued cell velocity equals the VoF velocity), so cut
          // cells are no longer processed twice and the min-reduction is unchanged.
          VoFIterator& vofit = (*m_amr->getMultiCutVofIterator(m_realm, m_phase)[lvl])[din];

          // Regular grid data.
          const BaseFab<Real>& veloReg    = velo.getSingleValuedFAB();
          const BaseFab<Real>& notCovered = (*m_amr->getNotCoveredCells(m_realm, m_phase)[lvl])[din]
                                              .getSingleValuedFAB();

          // Compute dt = dx/(|vx|+|vy|+|vz|) and check if it's smaller than the smallest so far.
          auto regularKernel = [&](const IntVect& iv) -> void {
            Real velMax = 0.0;
            if (notCovered(iv, m_comp) > 0.0) {
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
          BoxLoops::loop<D_DECL(1, 1, 1)>(cellBox, regularKernel);
          BoxLoops::loop(vofit, irregularKernel);
        }
      }
    }
  }

  return ParallelOps::min(minDt);
}

void
CdrCTU::parseSlopeLimiter()
{
  CH_TIME("CdrCTU::parseSlopeLimiter()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseSlopeLimiter()" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string slopeLimiter;
  pp.get("use_ctu", m_useCTU);
  pp.get("slope_limiter", slopeLimiter);
  pp.get("eb_inflow_state", m_ebInflowState);
  pp.get("flux_correction", m_fluxCorrection);

  if (slopeLimiter == "none") {
    m_limiter = Limiter::None;
  }
  else if (slopeLimiter == "minmod") {
    m_limiter = Limiter::MinMod;
  }
  else if (slopeLimiter == "superbee") {
    m_limiter = Limiter::Superbee;
  }
  else if (slopeLimiter == "mc") {
    m_limiter = Limiter::MonotonizedCentral;
  }
  else {
    MayDay::Error("CdrCTU::parseSlopeLimiter -- unknown limiter requested");
  }
}

void
CdrCTU::advectToFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_cellPhi, const Real a_extrapDt, const Real a_dt)
{
  CH_TIME("CdrCTU::advectToFaces(EBAMRFluxData, EBAMRCellData, Real, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::advectToFaces(EBAMRFluxData, EBAMRCellData, Real, Real)" << endl;
  }

  const int numberOfGhostCells = m_amr->getNumberOfGhostCells();

  CH_assert(a_facePhi[0]->nComp() == 1);
  CH_assert(a_cellPhi[0]->nComp() == 1);
  CH_assert(numberOfGhostCells >= 2);

  DataOps::setValue(a_facePhi, 0.0);

  // Ghost cells need to be interpolated. We make a copy of a_cellPhi which we use for that. This requires
  EBAMRCellData phi;
  m_amr->allocate(phi, m_realm, m_phase, m_nComp);
  DataOps::copy(phi, a_cellPhi);

  m_amr->conservativeAverage(phi, m_realm, m_phase);
  m_amr->interpGhostPwl(phi, m_realm, m_phase);

  // The EB inflow state on the faces of the grids needs the EB flux in the ghost cells.
  if (m_ebInflowState) {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      m_ebFlux[lvl]->exchange();
    }
  }

  // With flux correction, the first-order states (zero slopes) are kept for the correction.
  const bool useFluxCorrection = m_fluxCorrection && a_dt > 0.0 && m_limiter != Limiter::None;

  EBAMRFluxData facePhiLow;

  if (useFluxCorrection) {
    m_amr->allocate(facePhiLow, m_realm, m_phase, m_nComp);
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const ProblemDomain&     domain = m_amr->getDomains()[lvl];
    const DataIterator&      dit    = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      EBFluxFAB&       facePhi = (*a_facePhi[lvl])[din];
      const EBCellFAB& cellPhi = (*phi[lvl])[din];
      const EBCellFAB& cellVel = (*m_cellVelocity[lvl])[din];
      const EBFluxFAB& faceVel = (*m_faceVelocity[lvl])[din];

      const Box      cellBox = dbl[din];
      const EBISBox& ebisbox = ebisl[din];

      // Limit slopes and solve Riemann problem (which yields the upwind state at the face). Note that we need one ghost
      // cell for the slopes because in order to extrapolate to the left/right sides of a face, we need the centered
      // slope on both sides for the upwind. So, normalSlopes is bigger than cellBox (by one). Since the limited slope
      // is computed using the left/right slopes, we end up needing two grid cells.
      Box grownBox = cellBox;
      grownBox.grow(1);
      EBCellFAB normalSlopes(ebisbox, grownBox, SpaceDim);
      normalSlopes.setVal(0.0);

      // Compute normal slopes.
      if (m_limiter != Limiter::None) {
        this->computeNormalSlopes(normalSlopes, cellPhi, cellBox, domain, lvl, din);
      }

      this->upwind(facePhi, normalSlopes, cellPhi, cellVel, faceVel, domain, cellBox, lvl, din, a_extrapDt);

      if (m_ebInflowState && a_dt > 0.0) {
        for (int dir = 0; dir < SpaceDim; dir++) {
          this->ebInflowState(facePhi[dir], cellVel, faceVel[dir], cellBox, dir, lvl, din, a_dt);
        }
      }

      if (useFluxCorrection) {
        EBFluxFAB& facePhiFirst = (*facePhiLow[lvl])[din];
        EBCellFAB  zeroSlopes(ebisbox, grownBox, SpaceDim);

        zeroSlopes.setVal(0.0);

        this->upwind(facePhiFirst, zeroSlopes, cellPhi, cellVel, faceVel, domain, cellBox, lvl, din, a_extrapDt);

        if (m_ebInflowState) {
          for (int dir = 0; dir < SpaceDim; dir++) {
            this->ebInflowState(facePhiFirst[dir], cellVel, faceVel[dir], cellBox, dir, lvl, din, a_dt);
          }
        }
      }
    }
  }

  if (useFluxCorrection) {
    this->fluxCorrection(a_facePhi, facePhiLow, phi, a_dt);
  }
}

void
CdrCTU::computeNormalSlopes(EBCellFAB&           a_normalSlopes,
                            const EBCellFAB&     a_cellPhi,
                            const Box&           a_cellBox,
                            const ProblemDomain& a_domain,
                            const int            a_level,
                            const DataIndex&     a_dit)
{
  CH_TIME("CdrCTU::computeNormalSlopes(EBCellFAB, EBCellFAB, Box, ProblemDomain, int, DataIndex)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeNormalSlopes(EBCellFAB, EBCellFAB, Box, ProblemDomain, int, DataIndex)" << endl;
  }

  CH_assert(a_normalSlopes.nComp() == SpaceDim);
  CH_assert(a_cellPhi.nComp() == 1);

  const Box&     domainBox = a_domain.domainBox();
  const EBISBox& ebisbox   = a_cellPhi.getEBISBox();

  // Compute slopes in regular cells.
  BaseFab<Real>&       slopesReg = a_normalSlopes.getSingleValuedFAB();
  const BaseFab<Real>& phiReg    = a_cellPhi.getSingleValuedFAB();

  // Compute slopes in regular grid cells. Note that we need to grow the input box by one ghost cell in direction 'dir'
  // because we need the centered slopes on both sides of the faces that we extrapolate to.
  for (int dir = 0; dir < SpaceDim; dir++) {

    // Grown box
    const Box grownBox = grow(a_cellBox, 1) & domainBox;

    // Computation box for regular kernel; we need slopes outside
    // the grid patch so grow the box by one, but don't include
    // boundary cells.
    Box interiorCells = grow(a_cellBox, 1);
    interiorCells &= grow(domainBox, -1);

    // Boundary cells on the low/high sides.
    Box bndryLo = adjCellLo(domainBox, dir, -1);
    Box bndryHi = adjCellHi(domainBox, dir, -1);

    bndryLo &= grownBox;
    bndryHi &= grownBox;

    CH_assert(a_domain.contains(interiorCells));
    CH_assert(a_domain.contains(bndryLo));
    CH_assert(a_domain.contains(bndryHi));

    // Irregular kernel domain -- also does boundary cells if they are cut.
    VoFIterator& vofit = (*m_grownVofIter[a_level])[a_dit];

    const IntVect shift = BASISV(dir);

    // Set up the regular slope kernels. Add more slopes if you need them.
    std::function<void(const IntVect&)> regularKernel;
    switch (m_limiter) {
    case Limiter::None: {
      regularKernel = [&](const IntVect& iv) -> void {
        CH_assert(a_domain.contains(iv));
        slopesReg(iv, dir) = 0.0;
      };

      break;
    }
    case Limiter::MinMod: {
      regularKernel = [&](const IntVect& iv) -> void {
        CH_assert(a_domain.contains(iv + shift));
        CH_assert(a_domain.contains(iv - shift));

        const Real dwl = phiReg(iv, m_comp) - phiReg(iv - shift, m_comp);
        const Real dwr = phiReg(iv + shift, m_comp) - phiReg(iv, m_comp);

        if (dwl * dwr > 0.0) {
          slopesReg(iv, dir) = ChomboDischarge::CdrCTU::minmod(dwl, dwr);
        }
        else {
          slopesReg(iv, dir) = 0.0;
        }
      };

      break;
    }
    case Limiter::Superbee: {
      regularKernel = [&](const IntVect& iv) -> void {
        CH_assert(a_domain.contains(iv + shift));
        CH_assert(a_domain.contains(iv - shift));

        const Real dwl = phiReg(iv, m_comp) - phiReg(iv - shift, m_comp);
        const Real dwr = phiReg(iv + shift, m_comp) - phiReg(iv, m_comp);

        if (dwl * dwr > 0.0) {
          slopesReg(iv, dir) = superbee(dwl, dwr);
        }
        else {
          slopesReg(iv, dir) = 0.0;
        }
      };

      break;
    }
    case Limiter::MonotonizedCentral: {
      regularKernel = [&](const IntVect& iv) -> void {
        const Real dwl = phiReg(iv, m_comp) - phiReg(iv - shift, m_comp);
        const Real dwr = phiReg(iv + shift, m_comp) - phiReg(iv, m_comp);

        if (dwl * dwr > 0.0) {
          slopesReg(iv, dir) = ChomboDischarge::CdrCTU::monotonizedCentral(dwl, dwr);
        }
        else {
          slopesReg(iv, dir) = 0.0;
        }
      };

      break;
    }
    default: {
      MayDay::Error("CdrCTU::computeNormalSlopes -- logic bust");

      break;
    }
    }

    // Kernel for cells abutting the boundaries.
    auto boundaryKernelLo = [&](const IntVect& iv) -> void {
      CH_assert(a_domain.contains(iv));
      CH_assert(a_domain.contains(iv + shift));

      slopesReg(iv, dir) = phiReg(iv + shift, m_comp) - phiReg(iv, m_comp);
    };

    auto boundaryKernelHi = [&](const IntVect& iv) -> void {
      CH_assert(a_domain.contains(iv));
      CH_assert(a_domain.contains(iv - shift));

      slopesReg(iv, dir) = phiReg(iv, m_comp) - phiReg(iv - shift, m_comp);
    };

    // Cut-cell kernel -- also does cut-cells on the boundary.
    auto irregularKernel = [&](const VolIndex& vof) -> void {
      const IntVect iv = vof.gridIndex();

      const bool onLoSide = (iv[dir] == domainBox.smallEnd(dir));
      const bool onHiSide = (iv[dir] == domainBox.bigEnd(dir));

      const bool hasFacesLeft = (ebisbox.numFaces(vof, dir, Side::Lo) == 1) && !onLoSide;
      const bool hasFacesRigh = (ebisbox.numFaces(vof, dir, Side::Hi) == 1) && !onHiSide;

      VolIndex vofLeft;
      VolIndex vofRigh;

      Real dwl = 0.0;
      Real dwr = 0.0;

      Real phiLeft = 0.0;
      Real phiRigh = 0.0;

      // Compute left and right slope
      if (hasFacesLeft) {
        Vector<FaceIndex> facesLeft = ebisbox.getFaces(vof, dir, Side::Lo);
        vofLeft                     = facesLeft[0].getVoF(Side::Lo);
        phiLeft                     = a_cellPhi(vofLeft, m_comp);
        dwl                         = a_cellPhi(vof, m_comp) - phiLeft;
      }
      if (hasFacesRigh) {
        Vector<FaceIndex> facesRigh = ebisbox.getFaces(vof, dir, Side::Hi);
        vofRigh                     = facesRigh[0].getVoF(Side::Hi);
        phiRigh                     = a_cellPhi(vofRigh, m_comp);
        dwr                         = phiRigh - a_cellPhi(vof, m_comp);
      }

      if (!hasFacesLeft && hasFacesRigh) {
        dwl = dwr;
      }
      else if (hasFacesLeft && !hasFacesRigh) {
        dwr = dwl;
      }

      // Where the flow leaves the EB, the neighbour on the side that the EB closes (no face, and not the domain
      // boundary) is the EB inflow state, one cell away, instead of a copy of the difference on the other side.
      const bool closedLeft = (ebisbox.numFaces(vof, dir, Side::Lo) == 0) && !onLoSide;
      const bool closedRigh = (ebisbox.numFaces(vof, dir, Side::Hi) == 0) && !onHiSide;

      if (m_ebInflowState && ((closedLeft && hasFacesRigh) || (closedRigh && hasFacesLeft))) {
        Real phiEB = 0.0;

        if (this->getEbInflowState(phiEB, vof, ebisbox, a_level, a_dit)) {
          if (closedLeft) {
            dwl = a_cellPhi(vof, m_comp) - phiEB;
          }
          else {
            dwr = phiEB - a_cellPhi(vof, m_comp);
          }
        }
      }

      // Limit the slopes.
      switch (m_limiter) {
      case Limiter::None: {
        a_normalSlopes(vof, dir) = 0.0;

        break;
      }
      case Limiter::MinMod: {
        a_normalSlopes(vof, dir) = ChomboDischarge::CdrCTU::minmod(dwl, dwr);

        break;
      }
      case Limiter::Superbee: {
        a_normalSlopes(vof, dir) = superbee(dwl, dwr);

        break;
      }
      case Limiter::MonotonizedCentral: {
        a_normalSlopes(vof, dir) = ChomboDischarge::CdrCTU::monotonizedCentral(dwl, dwr);

        break;
      }
      default: {
        MayDay::Error("CD_CdrCTU::computeNormalSlopes -- logic bust 2");
      }
      }
    };

    // Apply the kernels. Beware of corrected slopes near the boundaries.
    BoxLoops::loop<D_DECL(1, 1, 1)>(interiorCells, regularKernel);
    BoxLoops::loop<D_DECL(1, 1, 1)>(bndryLo, boundaryKernelLo);
    BoxLoops::loop<D_DECL(1, 1, 1)>(bndryHi, boundaryKernelHi);
    BoxLoops::loop(vofit, irregularKernel);
  }
}

void
CdrCTU::upwind(EBFluxFAB&           a_facePhi,
               const EBCellFAB&     a_normalSlopes,
               const EBCellFAB&     a_cellPhi,
               const EBCellFAB&     a_cellVel,
               const EBFluxFAB&     a_faceVel,
               const ProblemDomain& a_domain,
               const Box&           a_cellBox,
               const int&           a_level,
               const DataIndex& /*a_dit*/,
               const Real& a_dt)
{
  CH_TIME("CdrCTU::upwind(EBFluxFAB, EBCellFABx3, EBFluxFAB, ProblemDomain, Box, int, DataIndex, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::upwind(EBFluxFAB, EBCellFABx3, EBFluxFAB, ProblemDomain, Box, int, DataIndex, Real)" << endl;
  }

  CH_assert(a_facePhi.nComp() == 1);
  CH_assert(a_normalSlopes.nComp() == SpaceDim);
  CH_assert(a_cellPhi.nComp() == 1);
  CH_assert(a_cellVel.nComp() == SpaceDim);
  CH_assert(a_faceVel.nComp() == 1);

  // TLDR: We want to compute the states at cell-centers, and possible at the half time step. E.g. we
  //       want phi_[(i+1/2),j,k]^(k+1/2). The normal slopes (slopes that are normal to the face normal)
  //       come into this routine through a_normalSlopes, and the transverse slopes are added in this method.

  const Real dt  = m_useCTU ? a_dt : 0.0;
  const Real dx  = m_amr->getDx()[a_level];
  const Real dtx = dt / dx;

  const Box& domainBox = a_domain.domainBox();

  for (int dir = 0; dir < SpaceDim; dir++) {
    EBFaceFAB&       facePhi = a_facePhi[dir];
    const EBFaceFAB& faceVel = a_faceVel[dir];
    const EBISBox&   ebisbox = a_cellPhi.getEBISBox();
    const EBGraph&   ebgraph = ebisbox.getEBGraph();

    // Single-valued data.
    BaseFab<Real>&       regFacePhi = a_facePhi[dir].getSingleValuedFAB();
    const BaseFab<Real>& regSlopes  = a_normalSlopes.getSingleValuedFAB();
    const BaseFab<Real>& regStates  = a_cellPhi.getSingleValuedFAB();
    const BaseFab<Real>& regCellVel = a_cellVel.getSingleValuedFAB();
    const BaseFab<Real>& regFaceVel = a_faceVel[dir].getSingleValuedFAB();

    // Iteration space for the kernels. When upwinding we want to set phi on the faces, so the iteration space is
    // defined by the faces of this box (in direction dir).
    const Vector<FaceIndex> irregFaces = ebgraph.getIrregFaces(a_cellBox, dir);

    // Interior faces.
    Box interiorFaces = grow(a_cellBox, 1);
    interiorFaces &= a_domain;
    interiorFaces.grow(dir, -1);
    interiorFaces.surroundingNodes(dir);

    // Boundary faces on lo/high side.
    Box bndryFacesLo = adjCellLo(domainBox, dir, -1);
    Box bndryFacesHi = adjCellHi(domainBox, dir, -1);

    bndryFacesLo &= a_cellBox;
    bndryFacesHi &= a_cellBox;

    // In 3D, the corner term of the corner transport upwind method,
    //
    //    (dt/dx)^2 * |v_t| * |v_s| * (phi - phi_t - phi_s + phi_ts)/3,
    //
    // with t and s the transverse directions and phi_t, phi_s, phi_ts the upwind neighbours in t, in s, and across the
    // corner. Summed over the three face directions, it completes the update to the product of one-dimensional upwind
    // updates. Without it the update is unstable for Courant numbers above 1/2 in all three directions.
    auto regularCornerTerm = [&](const IntVect& a_cell) -> Real {
      Real term = 0.0;

#if CH_SPACEDIM == 3
      const int  t  = (dir + 1) % 3;
      const int  s  = (dir + 2) % 3;
      const Real vt = regCellVel(a_cell, t);
      const Real vs = regCellVel(a_cell, s);

      if (vt != 0.0 && vs != 0.0) {
        const IntVect et = (vt > 0.0) ? -BASISV(t) : BASISV(t);
        const IntVect es = (vs > 0.0) ? -BASISV(s) : BASISV(s);

        // Covered cells have no state, so the term is left out next to them.
        if (ebisbox.isCovered(a_cell + et) || ebisbox.isCovered(a_cell + es) || ebisbox.isCovered(a_cell + et + es)) {
          return term;
        }

        term = dtx * dtx * std::abs(vt) * std::abs(vs) *
               (regStates(a_cell, m_comp) - regStates(a_cell + et, m_comp) - regStates(a_cell + es, m_comp) +
                regStates(a_cell + et + es, m_comp)) /
               3.0;
      }
#endif

      return term;
    };

    // The upwind neighbour of a cut cell in a direction, if it is reached through a single face inside the domain.
    auto upwindNeighbour =
      [&](VolIndex& a_neighbour, const VolIndex& a_vof, const int a_dir, const Real a_vel) -> bool {
      const Side::LoHiSide side = (a_vel > 0.0) ? Side::Lo : Side::Hi;
      const IntVect        iv   = a_vof.gridIndex() + ((a_vel > 0.0) ? -BASISV(a_dir) : BASISV(a_dir));

      if (!a_domain.contains(iv)) {
        return false;
      }

      const Vector<FaceIndex> faces = ebisbox.getFaces(a_vof, a_dir, side);

      if (faces.size() != 1) {
        return false;
      }

      a_neighbour = faces[0].getVoF(side);

      return true;
    };

    // The corner term in cut cells. It is left out where the upwind neighbours cannot be reached.
    auto irregularCornerTerm = [&](const VolIndex& a_vof) -> Real {
      Real term = 0.0;

#if CH_SPACEDIM == 3
      const int  t  = (dir + 1) % 3;
      const int  s  = (dir + 2) % 3;
      const Real vt = a_cellVel(a_vof, t);
      const Real vs = a_cellVel(a_vof, s);

      VolIndex vofT;
      VolIndex vofS;
      VolIndex vofTS;

      if (vt != 0.0 && vs != 0.0 && upwindNeighbour(vofT, a_vof, t, vt) && upwindNeighbour(vofS, a_vof, s, vs) &&
          upwindNeighbour(vofTS, vofT, s, vs)) {
        term = dtx * dtx * std::abs(vt) * std::abs(vs) *
               (a_cellPhi(a_vof, m_comp) - a_cellPhi(vofT, m_comp) - a_cellPhi(vofS, m_comp) +
                a_cellPhi(vofTS, m_comp)) /
               3.0;
      }
#endif

      return term;
    };

    // Regular upwind kernel. This also includes the transverse terms; also not that we call this on the
    // face-centered box so for a grid face with index iv, the cell on the left side of the face has index
    // iv - BASISV(dir)
    auto regularKernel = [&](const IntVect& iv) -> void {
      // Cells that are left/right of the current face.
      const IntVect cellLeft = iv - BASISV(dir);
      const IntVect cellRigh = iv;

      // Normal extrapolation.
      Real primLeft = regStates(cellLeft, m_comp) +
                      0.5 * std::min(1.0, 1.0 - regCellVel(cellLeft, dir) * dtx) * regSlopes(cellLeft, dir);
      Real primRigh = regStates(cellRigh, m_comp) -
                      0.5 * std::min(1.0, 1.0 + regCellVel(cellRigh, dir) * dtx) * regSlopes(cellRigh, dir);

      // Compute the transverse (CTU) terms.
      for (int transverseDir = 0; transverseDir < SpaceDim; transverseDir++) {
        if (transverseDir != dir) {
          Real slopeLeft = 0.0;
          Real slopeRigh = 0.0;

          // Transverse term in cell to the left.
          if (regCellVel(cellLeft, transverseDir) < 0.0) {
            slopeLeft = regStates(cellLeft + BASISV(transverseDir), m_comp) - regStates(cellLeft, m_comp);
          }
          else if (regCellVel(cellLeft, transverseDir) > 0.0) {
            slopeLeft = regStates(cellLeft, m_comp) - regStates(cellLeft - BASISV(transverseDir), m_comp);
          }

          // Transverse term in cell to the right.
          if (regCellVel(cellRigh, transverseDir) < 0.0) {
            slopeRigh = regStates(cellRigh + BASISV(transverseDir), m_comp) - regStates(cellRigh, m_comp);
          }
          else if (regCellVel(cellRigh, transverseDir) > 0.0) {
            slopeRigh = regStates(cellRigh, m_comp) - regStates(cellRigh - BASISV(transverseDir), m_comp);
          }

          primLeft -= 0.5 * dtx * regCellVel(cellLeft, transverseDir) * slopeLeft;
          primRigh -= 0.5 * dtx * regCellVel(cellRigh, transverseDir) * slopeRigh;
        }
      }

      primLeft += regularCornerTerm(cellLeft);
      primRigh += regularCornerTerm(cellRigh);

      // Solve the Riemann problem.
      Real&       facePhi = regFacePhi(iv, m_comp);
      const Real& faceVel = regFaceVel(iv, m_comp);

      if (faceVel > 0.0) {
        facePhi = primLeft;
      }
      else if (faceVel < 0.0) {
        facePhi = primRigh;
      }
      else {
        facePhi = 0.0;
      }
    };

    // Transverse (CTU) terms of a cell, as in the regular kernel, added to an extrapolated state.
    auto addTransverseTerms = [&](Real& a_prim, const IntVect& a_cell) -> void {
      for (int transverseDir = 0; transverseDir < SpaceDim; transverseDir++) {
        if (transverseDir != dir) {
          Real slope = 0.0;

          if (regCellVel(a_cell, transverseDir) < 0.0) {
            slope = regStates(a_cell + BASISV(transverseDir), m_comp) - regStates(a_cell, m_comp);
          }
          else if (regCellVel(a_cell, transverseDir) > 0.0) {
            slope = regStates(a_cell, m_comp) - regStates(a_cell - BASISV(transverseDir), m_comp);
          }

          a_prim -= 0.5 * dtx * regCellVel(a_cell, transverseDir) * slope;
        }
      }

      a_prim += regularCornerTerm(a_cell);
    };

    // Kernel for cells abutting the boundaries on the low end. Note that the input box
    // to this kernel is a cell-centered box, so the input IntVect is the
    // both the face index and the cell index to the right of the face.
    auto boundaryKernelLo = [&](const IntVect& iv) -> void {
      // Normal extrapolation from the cell to the right of the face, and the transverse terms.
      Real primLeft = 0.0;
      Real primRigh = regStates(iv, m_comp) - 0.5 * std::min(1.0, 1.0 + regCellVel(iv, dir) * dtx) * regSlopes(iv, dir);

      addTransverseTerms(primRigh, iv);

      // Solve the Riemann problem.
      Real&       facePhi = regFacePhi(iv, m_comp);
      const Real& faceVel = regFaceVel(iv, m_comp);

      if (faceVel > 0.0) {
        facePhi = primLeft;
      }
      else if (faceVel < 0.0) {
        facePhi = primRigh;
      }
      else {
        facePhi = 0.0;
      }
    };

    // Kernel for cells abutting the boundaries on the low end. Note that the input box
    // to this kernel is a cell-centered box, so the boundary face index is the input
    // IntVect + BASISV(dir)
    // both the face index and the cell index to the right of the face.
    auto boundaryKernelHi = [&](const IntVect& iv) -> void {
      // Normal extrapolation from the cell to the left of the face, and the transverse terms.
      Real primLeft = regStates(iv, m_comp) + 0.5 * std::min(1.0, 1.0 - regCellVel(iv, dir) * dtx) * regSlopes(iv, dir);
      Real primRigh = 0.0;

      addTransverseTerms(primLeft, iv);

      // Solve the Riemann problem.
      Real&       facePhi = regFacePhi(iv + BASISV(dir), m_comp);
      const Real& faceVel = regFaceVel(iv + BASISV(dir), m_comp);

      if (faceVel > 0.0) {
        facePhi = primLeft;
      }
      else if (faceVel < 0.0) {
        facePhi = primRigh;
      }
      else {
        facePhi = 0.0;
      }
    };

    // Cut-cell kernel. Same as the above but we need to explicitly ensure that
    // we are getting the correct left/right vofs (cells might be multi-valued).
    auto irregularKernel = [&](const FaceIndex& face) -> void {
      const VolIndex& vofLeft = face.getVoF(Side::Lo);
      const VolIndex& vofRigh = face.getVoF(Side::Hi);

      if (!face.isBoundary()) {
        CH_assert(a_domain.contains(vofLeft.gridIndex()));
        CH_assert(a_domain.contains(vofRigh.gridIndex()));

        Real primLeft = a_cellPhi(vofLeft, m_comp) +
                        0.5 * std::min(1.0, 1.0 - a_cellVel(vofLeft, dir) * dtx) * a_normalSlopes(vofLeft, dir);
        Real primRigh = a_cellPhi(vofRigh, m_comp) -
                        0.5 * std::min(1.0, 1.0 + a_cellVel(vofRigh, dir) * dtx) * a_normalSlopes(vofRigh, dir);

        // Compute the transverse (CTU) terms.
        if (m_useCTU) {
          for (int transverseDir = 0; transverseDir < SpaceDim; transverseDir++) {

            if (transverseDir != dir) {
              Real slopeLeft = 0.0;
              Real slopeRigh = 0.0;

              // Availability of cells used for transverse slopes. We turn off the transverse
              // terms if the required cells are not available.
              const IntVect ivLeftUp   = vofLeft.gridIndex() + BASISV(transverseDir);
              const IntVect ivLeftDown = vofLeft.gridIndex() - BASISV(transverseDir);
              const IntVect ivRighUp   = vofRigh.gridIndex() + BASISV(transverseDir);
              const IntVect ivRighDown = vofRigh.gridIndex() - BASISV(transverseDir);

              // Transverse term in cell to the left.
              if (a_cellVel(vofLeft, transverseDir) < 0.0 && a_domain.contains(ivLeftUp)) {
                const Vector<FaceIndex>& facesHi = ebisbox.getFaces(vofLeft, transverseDir, Side::Hi);

                const int nFaces = static_cast<int>(facesHi.size());

                if (nFaces > 0) {
                  Real phiUp = 0.0;
                  for (int iface = 0; iface < nFaces; iface++) {
                    phiUp += a_cellPhi(facesHi[iface].getVoF(Side::Hi), m_comp);
                  }

                  phiUp /= nFaces;

                  slopeLeft = phiUp - a_cellPhi(vofLeft, m_comp);
                }
              }
              else if (a_cellVel(vofLeft, transverseDir) > 0.0 && a_domain.contains(ivLeftDown)) {
                const Vector<FaceIndex>& facesLo = ebisbox.getFaces(vofLeft, transverseDir, Side::Lo);

                const int nFaces = static_cast<int>(facesLo.size());

                if (nFaces > 0) {
                  Real phiDown = 0.0;
                  for (int iface = 0; iface < nFaces; iface++) {
                    phiDown += a_cellPhi(facesLo[iface].getVoF(Side::Lo), m_comp);
                  }

                  phiDown /= nFaces;

                  slopeLeft = a_cellPhi(vofLeft, m_comp) - phiDown;
                }
              }

              // Transverse term in cell to the right
              if (a_cellVel(vofRigh, transverseDir) < 0.0 && a_domain.contains(ivRighUp)) {
                const Vector<FaceIndex>& facesHi = ebisbox.getFaces(vofRigh, transverseDir, Side::Hi);

                const int nFaces = static_cast<int>(facesHi.size());

                if (nFaces > 0) {
                  Real phiUp = 0.0;
                  for (int iface = 0; iface < nFaces; iface++) {
                    phiUp += a_cellPhi(facesHi[iface].getVoF(Side::Hi), m_comp);
                  }
                  phiUp /= nFaces;

                  slopeRigh = phiUp - a_cellPhi(vofRigh, m_comp);
                }
              }
              else if (a_cellVel(vofRigh, transverseDir) > 0.0 && a_domain.contains(ivRighDown)) {
                const Vector<FaceIndex>& facesLo = ebisbox.getFaces(vofRigh, transverseDir, Side::Lo);

                const int nFaces = static_cast<int>(facesLo.size());

                if (nFaces > 0) {
                  Real phiDown = 0.0;
                  for (int iface = 0; iface < nFaces; iface++) {
                    phiDown += a_cellPhi(facesLo[iface].getVoF(Side::Lo), m_comp);
                  }

                  phiDown /= nFaces;

                  slopeRigh = a_cellPhi(vofRigh, m_comp) - phiDown;
                }
              }

              primLeft -= 0.5 * dtx * a_cellVel(vofLeft, transverseDir) * slopeLeft;
              primRigh -= 0.5 * dtx * a_cellVel(vofRigh, transverseDir) * slopeRigh;
            }
          }

          primLeft += irregularCornerTerm(vofLeft);
          primRigh += irregularCornerTerm(vofRigh);
        }

        // Solve the Riemann problem.
        if (faceVel(face, m_comp) > 0.0) {
          facePhi(face, m_comp) = primLeft;
        }
        else if (faceVel(face, m_comp) < 0.0) {
          facePhi(face, m_comp) = primRigh;
        }
        else {
          facePhi(face, m_comp) = 0.0;
        }
      }
      else {
        // On the low/high boundary we don't have a characteristic coming into the domain so just set
        // the flux to whatever. If the domain face is multi-valued then set the outgoing flux equal to
        // the single-valued face flux.
        if (a_domain.contains(vofRigh.gridIndex())) {
          // Lo boundary face.
          facePhi(face, m_comp) = regFacePhi(vofRigh.gridIndex());
        }
        if (a_domain.contains(vofLeft.gridIndex())) {
          facePhi(face, m_comp) = regFacePhi(vofLeft.gridIndex() + BASISV(dir));
        }
      }
    };

    // Launch the kernels.
    BoxLoops::loop<D_DECL(1, 1, 1)>(interiorFaces, regularKernel);
    BoxLoops::loop<D_DECL(1, 1, 1)>(bndryFacesLo, boundaryKernelLo);
    BoxLoops::loop<D_DECL(1, 1, 1)>(bndryFacesHi, boundaryKernelHi);
    BoxLoops::loop(irregFaces, irregularKernel);
  }
}

void
CdrCTU::fluxCorrection(EBAMRFluxData&       a_facePhi,
                       const EBAMRFluxData& a_facePhiLow,
                       const EBAMRCellData& a_cellPhi,
                       const Real           a_dt) noexcept
{
  CH_TIME("CdrCTU::fluxCorrection");
  if (m_verbosity > 5) {
    pout() << m_name + "::fluxCorrection" << endl;
  }

  // Fraction R of the correction that may take mass from each cell. It is 1 where nothing is limited, also in ghost
  // cells that no box owns.
  EBAMRCellData fraction;
  m_amr->allocate(fraction, m_realm, m_phase, 1, 1);
  DataOps::setValue(fraction, 1.0);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const DataIterator&      dit   = dbl.dataIterator();
    const Real               dtx   = a_dt / m_amr->getDx()[lvl];

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din     = dit[mybox];
      const Box        box     = dbl[din];
      const EBISBox&   ebisbox = ebisl[din];

      EBCellFAB&             frac    = (*fraction[lvl])[din];
      const EBCellFAB&       phi     = (*a_cellPhi[lvl])[din];
      const EBFluxFAB&       phiHigh = (*a_facePhi[lvl])[din];
      const EBFluxFAB&       phiLow  = (*a_facePhiLow[lvl])[din];
      const EBFluxFAB&       faceVel = (*m_faceVelocity[lvl])[din];
      const BaseIVFAB<Real>& ebFlux  = (*m_ebFlux[lvl])[din];

      BaseFab<Real>&       fracReg = frac.getSingleValuedFAB();
      const BaseFab<Real>& phiReg  = phi.getSingleValuedFAB();

      // Regular cells: the first-order update, and the mass that the correction takes out of the cell. Cut cells are
      // done again below with their volume and area fractions and the EB flux.
      auto regularKernel = [&](const IntVect& iv) -> void {
        Real netLow  = 0.0;
        Real removed = 0.0;

        for (int dir = 0; dir < SpaceDim; dir++) {
          const BaseFab<Real>& fh = phiHigh[dir].getSingleValuedFAB();
          const BaseFab<Real>& fl = phiLow[dir].getSingleValuedFAB();
          const BaseFab<Real>& fv = faceVel[dir].getSingleValuedFAB();

          const IntVect hi = iv + BASISV(dir);

          netLow += fv(hi, m_comp) * fl(hi, m_comp) - fv(iv, m_comp) * fl(iv, m_comp);
          removed += std::max(0.0, fv(hi, m_comp) * (fh(hi, m_comp) - fl(hi, m_comp)));
          removed += std::max(0.0, -fv(iv, m_comp) * (fh(iv, m_comp) - fl(iv, m_comp)));
        }

        const Real massLow = std::max(0.0, phiReg(iv, m_comp) - dtx * netLow);

        removed *= dtx;

        fracReg(iv, m_comp) = (removed > massLow) ? massLow / removed : 1.0;
      };

      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[din];

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        Real netLow  = ebFlux(vof, m_comp) * ebisbox.bndryArea(vof);
        Real removed = 0.0;

        for (int dir = 0; dir < SpaceDim; dir++) {
          for (SideIterator sit; sit.ok(); ++sit) {
            const int               isign = sign(sit());
            const Vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, sit());

            for (int iface = 0; iface < faces.size(); iface++) {
              const FaceIndex& face = faces[iface];
              const Real       area = ebisbox.areaFrac(face);
              const Real       v    = faceVel[dir](face, m_comp);

              netLow += isign * area * v * phiLow[dir](face, m_comp);
              removed += std::max(0.0, isign * area * v * (phiHigh[dir](face, m_comp) - phiLow[dir](face, m_comp)));
            }
          }
        }

        const Real massLow = std::max(0.0, ebisbox.volFrac(vof) * phi(vof, m_comp) - dtx * netLow);

        removed *= dtx;

        frac(vof, m_comp) = (removed > massLow) ? massLow / removed : 1.0;
      };

      BoxLoops::loop(box, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }

    fraction[lvl]->exchange();

    // The corrected face states. The correction on a face takes mass from the cell on one side and is limited by the
    // fraction of that cell.
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din     = dit[mybox];
      const Box        box     = dbl[din];
      const EBISBox&   ebisbox = ebisl[din];

      const BaseFab<Real>& fracReg = (*fraction[lvl])[din].getSingleValuedFAB();

      for (int dir = 0; dir < SpaceDim; dir++) {
        BaseFab<Real>&       fh = (*a_facePhi[lvl])[din][dir].getSingleValuedFAB();
        const BaseFab<Real>& fl = (*a_facePhiLow[lvl])[din][dir].getSingleValuedFAB();
        const BaseFab<Real>& fv = (*m_faceVelocity[lvl])[din][dir].getSingleValuedFAB();

        auto faceKernel = [&](const IntVect& iv) -> void {
          const IntVect left = iv - BASISV(dir);

          if (ebisbox.isMultiValued(left) || ebisbox.isMultiValued(iv)) {
            return;
          }

          const Real correction = fh(iv, m_comp) - fl(iv, m_comp);
          const Real flux       = fv(iv, m_comp) * correction;

          Real c = 1.0;

          if (flux > 0.0) {
            c = fracReg(left, m_comp);
          }
          else if (flux < 0.0) {
            c = fracReg(iv, m_comp);
          }

          fh(iv, m_comp) = fl(iv, m_comp) + c * correction;
        };

        BoxLoops::loop(surroundingNodes(box, dir), faceKernel);
      }
    }
  }
}

RealVect
CdrCTU::getEbNormalIntoFluid(const VolIndex& a_vof, const EBISBox& a_ebisbox) const noexcept
{
  RealVect n = a_ebisbox.normal(a_vof);

  // Orient the normal from the EB centroid towards the volume centroid.
  if (n.dotProduct(a_ebisbox.centroid(a_vof) - a_ebisbox.bndryCentroid(a_vof)) < 0.0) {
    n = -n;
  }

  return n;
}

bool
CdrCTU::getEbInflowState(Real&            a_phiEB,
                         const VolIndex&  a_vof,
                         const EBISBox&   a_ebisbox,
                         const int        a_level,
                         const DataIndex& a_dit) const noexcept
{
  // The EB flux is only known in cells of this level's grids, not in ghost cells on refinement boundaries.
  if (a_ebisbox.bndryArea(a_vof) <= 0.0 || (*m_levelCells[a_level])[a_dit](a_vof.gridIndex(), 0) <= 0.0) {
    return false;
  }

  const RealVect   n       = this->getEbNormalIntoFluid(a_vof, a_ebisbox);
  const EBCellFAB& cellVel = (*m_cellVelocity[a_level])[a_dit];

  RealVect v;

  for (int dir = 0; dir < SpaceDim; dir++) {
    v[dir] = cellVel(a_vof, dir);
  }

  const Real vn = v.dotProduct(n);

  if (vn <= 0.0) {
    return false;
  }

  a_phiEB = -(*m_ebFlux[a_level])[a_dit](a_vof, m_comp) / vn;

  return true;
}

void
CdrCTU::ebInflowState(EBFaceFAB&       a_facePhi,
                      const EBCellFAB& a_cellVel,
                      const EBFaceFAB& a_faceVel,
                      const Box&       a_cellBox,
                      const int        a_dir,
                      const int        a_level,
                      const DataIndex& a_dit,
                      const Real       a_dt) noexcept
{
  CH_TIME("CdrCTU::ebInflowState");
  if (m_verbosity > 5) {
    pout() << m_name + "::ebInflowState" << endl;
  }

  const Real     dx      = m_amr->getDx()[a_level];
  const EBISBox& ebisbox = a_cellVel.getEBISBox();
  const Box      faceBox = surroundingNodes(a_cellBox, a_dir);

  // The cut cells of the box and its first layer of ghost cells, which are upwind of the faces of the box.
  VoFIterator& vofit = (*m_grownVofIter[a_level])[a_dit];

  auto kernel = [&](const VolIndex& vof) -> void {
    Real phiEB = 0.0;

    // Only where the flow leaves the EB.
    if (!this->getEbInflowState(phiEB, vof, ebisbox, a_level, a_dit)) {
      return;
    }

    const RealVect n = this->getEbNormalIntoFluid(vof, ebisbox);
    RealVect       v;

    for (int dir = 0; dir < SpaceDim; dir++) {
      v[dir] = a_cellVel(vof, dir);
    }

    const Real vn = v.dotProduct(n);

    for (SideIterator sit; sit.ok(); ++sit) {
      const Side::LoHiSide side  = sit();
      const int            isign = sign(side);

      const Vector<FaceIndex> faces = ebisbox.getFaces(vof, a_dir, side);

      for (int iface = 0; iface < faces.size(); iface++) {
        const FaceIndex& face = faces[iface];

        if (face.isBoundary() || !faceBox.contains(face.gridIndex(Side::Hi))) {
          continue;
        }

        // Outgoing faces only.
        if (isign * a_faceVel(face, m_comp) <= 0.0) {
          continue;
        }

        // Face centroid relative to the cell center, and its distance from the EB plane, in units of dx.
        RealVect xf = ebisbox.centroid(face);

        xf[a_dir] = 0.5 * isign;

        const Real s   = std::max(0.0, (xf - ebisbox.bndryCentroid(vof)).dotProduct(n));
        const Real tEB = s * dx / vn;
        const Real f   = std::min(1.0, tEB / a_dt);

        a_facePhi(face, m_comp) = f * a_facePhi(face, m_comp) + (1.0 - f) * phiEB;
      }
    }
  };

  BoxLoops::loop(vofit, kernel);
}

Real
CdrCTU::minmod(const Real& dwl, const Real& dwr) noexcept
{
  Real slope = 0.0;

  if (dwl * dwr > 0.0) {
    slope = std::abs(dwl) < std::abs(dwr) ? dwl : dwr;
  }

  return slope;
}

Real
CdrCTU::superbee(const Real& dwl, const Real& dwr) noexcept
{
  Real slope = 0.0;

  if (dwl * dwr > 0.0) {
    const Real s1 = ChomboDischarge::CdrCTU::minmod(dwl, 2 * dwr);
    const Real s2 = ChomboDischarge::CdrCTU::minmod(dwr, 2 * dwl);

    if (s1 * s2 > 0.0) {
      slope = std::abs(s1) > std::abs(s2) ? s1 : s2;
    }
  }

  return slope;
}

Real
CdrCTU::monotonizedCentral(const Real& dwl, const Real& dwr) noexcept
{
  Real slope = 0.0;

  if (dwl * dwr > 0.0) {
    const Real dwc = dwl + dwr;
    const Real sgn = Real((dwc > 0.0) - (dwc < 0.0));

    slope = sgn * std::min(0.5 * std::abs(dwc), 2.0 * std::min(std::abs(dwl), std::abs(dwr)));
  }

  return slope;
}

#include <CD_NamespaceFooter.H>
