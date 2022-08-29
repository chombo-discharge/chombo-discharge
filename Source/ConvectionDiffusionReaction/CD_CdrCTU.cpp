/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrCTU.cpp
  @brief  Implementation of CD_CdrCTU.H
  @author Robert Marskar
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
{
  CH_TIME("CdrCTU::CdrCTU()");

  // Class and object name
  m_className = "CdrCTU";
  m_name      = "CdrCTU";
  m_limiter   = Limiter::MonotonizedCentral;
  m_useCTU    = true;
}

CdrCTU::~CdrCTU() { CH_TIME("CdrCTU::~CdrCTU()"); }

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

    // TLDR: For advection, Bell, Collela, and Glaz says we must have dt <= dx/max(|vx|, |vy|, |vz|). See these three papers for details:
    //
    //       Colella, J. Comp. Phys. 87 (171-200), 1990
    //       Bell, Colella, Glaz, J. Comp. Phys 85 (257), 1989
    //       Minion, J. Comp. Phys 123 (435), 1996
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
  }

  return minDt;
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
CdrCTU::advectToFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_cellPhi, const Real a_dt)
{
  CH_TIME("CdrCTU::advectToFaces(EBAMRFluxData, EBAMRCellData, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::advectToFaces(EBAMRFluxData, EBAMRCellData, Real)" << endl;
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

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const ProblemDomain&     domain = m_amr->getDomains()[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      EBFluxFAB&       facePhi = (*a_facePhi[lvl])[dit()];
      const EBCellFAB& cellPhi = (*phi[lvl])[dit()];
      const EBCellFAB& cellVel = (*m_cellVelocity[lvl])[dit()];
      const EBFluxFAB& faceVel = (*m_faceVelocity[lvl])[dit()];

      const Box      cellBox = dbl[dit()];
      const EBISBox& ebisbox = ebisl[dit()];

      // Limit slopes and solve Riemann problem (which yields the upwind state at the face). Note that we need one ghost cell for
      // the slopes because in order to extrapolate to the left/right sides of a face, we need the centered slope on both
      // sides for the upwind. So, normalSlopes is bigger than cellBox (by one). Since the limited slope is computed using the
      // left/right slopes, we end up needing two grid cells.
      Box grownBox = cellBox;
      grownBox.grow(1);
      EBCellFAB normalSlopes(ebisbox, grownBox, SpaceDim);
      normalSlopes.setVal(0.0);

      // Compute normal slopes.
      if (m_limiter != Limiter::None) {
        this->computeNormalSlopes(normalSlopes, cellPhi, cellBox, domain, lvl, dit());
      }

      this->upwind(facePhi, normalSlopes, cellPhi, cellVel, faceVel, domain, cellBox, lvl, dit(), a_dt);
    }
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
  const EBGraph& ebgraph   = ebisbox.getEBGraph();

  // Compute slopes in regular cells.
  BaseFab<Real>&       slopesReg = a_normalSlopes.getSingleValuedFAB();
  const BaseFab<Real>& phiReg    = a_cellPhi.getSingleValuedFAB();

  // Compute slopes in regular grid cells. Note that we need to grow the input box by one ghost cell in direction 'dir' because we need the
  // centered slopes on both sides of the faces that we extrapolate to.
  for (int dir = 0; dir < SpaceDim; dir++) {

    // We don't want to compute slopes using data that lies outside the domain boundary -- that data can be bogus. So,
    // compute the strip of cells and check if these cells are in a_cellBox.
    const Box bndryLo = adjCellLo(domainBox, dir, 1); // Strip of cells on the low side in coordinate direction dir
    const Box bndryHi = adjCellHi(domainBox, dir, 1); // Strip of cells on the high side in coordinate direction dir

    CH_assert(domainBox.cellCentered());
    CH_assert(bndryLo.cellCentered());
    CH_assert(bndryHi.cellCentered());

    // We need the slopes outside the grid patch so grow the box by 1.
    const Box grownBox = grow(a_cellBox, 1) & domainBox;
    const Box loBox    = grownBox & bndryLo;
    const Box hiBox    = grownBox & bndryHi;

    // Computation box for regular kernel; we want to
    Box compBox = grow(a_cellBox, 1) & domainBox;

    compBox.shift(dir, -1);
    compBox &= domainBox;
    compBox.shift(dir, 1);

    compBox.shift(dir, 1);
    compBox &= domainBox;
    compBox.shift(dir, -1);

    // Irregular kernel domain -- we go through the same cells as in grownBox, but including only cut-cells.
    const IntVectSet irregIVS = ebisbox.getIrregIVS(grownBox);
    VoFIterator      vofit(irregIVS, ebgraph);

    const IntVect shift = BASISV(dir);

    // Set up the regular slope kernels. Add more slopes if you need them.
    std::function<void(const IntVect&)> regularKernel;
    switch (m_limiter) {
    case Limiter::None: {
      regularKernel = [&](const IntVect& iv) -> void { slopesReg(iv, dir) = 0.0; };

      break;
    }
    case Limiter::MinMod: {
      regularKernel = [&](const IntVect& iv) -> void {
        CH_assert(a_domain.contains(iv + shift));
        CH_assert(a_domain.contains(iv - shift));

        const Real dwl = phiReg(iv, m_comp) - phiReg(iv - shift, m_comp);
        const Real dwr = phiReg(iv + shift, m_comp) - phiReg(iv, m_comp);

        if (dwl * dwr > 0.0) {
          slopesReg(iv, dir) = this->minmod(dwl, dwr);
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
          slopesReg(iv, dir) = this->superbee(dwl, dwr);
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
          slopesReg(iv, dir) = this->monotonizedCentral(dwl, dwr);
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

    // Cut-cell kernel.
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

      // Limit the slopes.
      switch (m_limiter) {
      case Limiter::None: {
        a_normalSlopes(vof, dir) = 0.0;

        break;
      }
      case Limiter::MinMod: {
        a_normalSlopes(vof, dir) = this->minmod(dwl, dwr);

        break;
      }
      case Limiter::Superbee: {
        a_normalSlopes(vof, dir) = this->superbee(dwl, dwr);

        break;
      }
      case Limiter::MonotonizedCentral: {
        a_normalSlopes(vof, dir) = this->monotonizedCentral(dwl, dwr);

        break;
      }
      default: {
        MayDay::Error("CD_CdrCTU::computeNormalSlopes -- logic bust 2");
      }
      }
    };

    // Kernel for cells abutting the boundaries.
    auto boundaryKernelLo = [&](const IntVect& iv) -> void {
      slopesReg(iv, dir) = phiReg(iv + shift, m_comp) - phiReg(iv, m_comp);
    };
    auto boundaryKernelHi = [&](const IntVect& iv) -> void {
      slopesReg(iv, dir) = phiReg(iv, m_comp) - phiReg(iv - shift, m_comp);
    };

    // Apply the kernels. Beware of corrected slopes near the boundaries.
    BoxLoops::loop(compBox, regularKernel);

    if (!loBox.isEmpty()) {
      BoxLoops::loop(grownBox, boundaryKernelLo);
    }
    if (!hiBox.isEmpty()) {
      BoxLoops::loop(grownBox, boundaryKernelHi);
    }

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
               const DataIndex&     a_dit,
               const Real&          a_dt)
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

    // Iteration space for the kernels. When upwinding we want to set phi on the faces, so the iteration space is defined
    // by the faces of this box (in direction dir).
    //    const Box               faceBox    = surroundingNodes(a_cellBox, dir);
    const Vector<FaceIndex> irregFaces = ebgraph.getIrregFaces(a_cellBox, dir);

    // Clip boundary faces.
    Box faceBox = a_cellBox;
    faceBox.shift(dir, -1);
    faceBox &= a_domain;
    faceBox.shift(dir, 1);

    faceBox.shift(dir, 1);
    faceBox &= a_domain;
    faceBox.shift(dir, -1);

    faceBox.surroundingNodes(dir);

    // Regular upwind kernel. Recall that we call this on the faceBox, so the cell on the low side
    // has cell grid index iv - BASISV(dir)
    auto regularKernel = [&](const IntVect& iv) -> void {
      // Cells that are left/right of the current face.
      const IntVect cellLeft = iv - BASISV(dir);
      const IntVect cellRigh = iv;

      CH_assert(a_domain.contains(cellLeft));
      CH_assert(a_domain.contains(cellRigh));

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

    // Cut-cell kernel. Same as the above but we need to explicitly ensure that
    // we are getting the correct left/right vofs (cells might be multi-valued).
    auto irregularKernel = [&](const FaceIndex& face) -> void {
      if (!face.isBoundary()) {
        const VolIndex& vofLeft = face.getVoF(Side::Lo);
        const VolIndex& vofRigh = face.getVoF(Side::Hi);

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

                const int nFaces = facesHi.size();

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

                const int nFaces = facesLo.size();

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

                const int nFaces = facesHi.size();

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

                const int nFaces = facesLo.size();

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
    };

    // Launch the kernels.
    BoxLoops::loop(faceBox, regularKernel);
    BoxLoops::loop(irregFaces, irregularKernel);
  }
}

Real
CdrCTU::minmod(const Real& dwl, const Real& dwr) const noexcept
{
  Real slope = 0.0;

  if (dwl * dwr > 0.0) {
    slope = std::abs(dwl) < std::abs(dwr) ? dwl : dwr;
  }

  return slope;
}

Real
CdrCTU::superbee(const Real& dwl, const Real& dwr) const noexcept
{
  Real slope = 0.0;

  if (dwl * dwr > 0.0) {
    const Real s1 = this->minmod(dwl, 2 * dwr);
    const Real s2 = this->minmod(dwr, 2 * dwl);

    if (s1 * s2 > 0.0) {
      slope = std::abs(s1) > std::abs(s2) ? s1 : s2;
    }
  }

  return slope;
}

Real
CdrCTU::monotonizedCentral(const Real& dwl, const Real& dwr) const noexcept
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
