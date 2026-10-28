/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzEBStencils.cpp
 * @brief  Implementation of CD_PetscHelmholtzEBStencils.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>
#include <ParmParse.H>
#include <BoxIterator.H>

// Our includes
#include <CD_PetscHelmholtzEBStencils.H>
#include <CD_ParallelOps.H>
#include <CD_LeastSquares.H>
#include <CD_Location.H>
#include <CD_PolyUtils.H>
#include <CD_NamespaceHeader.H>

PetscHelmholtzEBStencils::PetscHelmholtzEBStencils() noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::PetscHelmholtzEBStencils");

  m_isDefined              = false;
  m_verbose                = false;
  m_order                  = -1;
  m_weight                 = -1;
  m_numCellsWithoutStencil = 0;

  ParmParse pp("PetscHelmholtzEBStencils");
  pp.query("verbose", m_verbose);
}

PetscHelmholtzEBStencils::~PetscHelmholtzEBStencils() noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::~PetscHelmholtzEBStencils");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBStencils::~PetscHelmholtzEBStencils" << endl;
  }
}

void
PetscHelmholtzEBStencils::setOrder(const int a_order) noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::setOrder");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBStencils::setOrder" << endl;
  }

  m_order = a_order;
}

void
PetscHelmholtzEBStencils::setWeight(const int a_weight) noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::setWeight");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBStencils::setWeight" << endl;
  }

  m_weight = a_weight;
}

void
PetscHelmholtzEBStencils::setInterface(const RefCountedPtr<BaseIF>& a_interface) noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::setInterface");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBStencils::setInterface" << endl;
  }

  m_interface = a_interface;
}

void
PetscHelmholtzEBStencils::clear() noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::clear");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBStencils::clear" << endl;
  }

  for (int iphase = 0; iphase < 2; iphase++) {
    m_derivStencils[iphase].resize(0);
    m_derivWeights[iphase].resize(0);
  }

  m_petscGrid              = RefCountedPtr<PetscGrid>();
  m_probLo                 = RealVect::Zero;
  m_numCellsWithoutStencil = 0;
  m_isDefined              = false;

  m_interpolators.resize(0);
  m_dx.resize(0);
}

void
PetscHelmholtzEBStencils::define(const RefCountedPtr<PetscGrid>&        a_petscGrid,
                                 const Vector<MFMultigridInterpolator>& a_interpolators,
                                 const RealVect&                        a_probLo,
                                 const Vector<Real>&                    a_dx) noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::define");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBStencils::define" << endl;
  }

  CH_assert(!a_petscGrid.isNull());
  CH_assert(a_petscGrid->isDefined());

  this->clear();

  const int numLevels = a_petscGrid->getMFLevelGrids().size();
  const int numPhases = a_petscGrid->getNumPhases();

  // Checked in every build: the stencils are indexed by these per level, and a bad order or weight gives no stencils.
  if (m_order < 1 || m_weight < 0) {
    MayDay::Abort("PetscHelmholtzEBStencils::define -- the order must be at least 1 and the weight non-negative");
  }
  if (a_dx.size() < numLevels || a_interpolators.size() < numLevels) {
    MayDay::Abort("PetscHelmholtzEBStencils::define -- dx or the interpolators do not cover every level");
  }

  // The neighbourhood radius is the order, and the interpolators fill only the ghost cells within their radius of a cut
  // cell. Checked in every build: a stencil reaching further would abort in the interpolator, deep in the assembly.
  for (int lvl = 1; lvl < numLevels; lvl++) {
    if (m_order > a_interpolators[lvl].getGhostCF()) {
      MayDay::Abort("PetscHelmholtzEBStencils::define -- the order exceeds the interpolation radius of the "
                    "coarse-fine interpolators (AmrMesh.mg_interp_radius)");
    }
  }

  m_petscGrid     = a_petscGrid;
  m_interpolators = a_interpolators;
  m_probLo        = a_probLo;
  m_dx            = a_dx;

  // Cut cells with an EB face but no accepted fit, which have no flux through their EB face.
  long long numWithoutStencil = 0;

  for (int iphase = 0; iphase < numPhases; iphase++) {
    m_derivStencils[iphase].resize(numLevels);
    m_derivWeights[iphase].resize(numLevels);

    for (int lvl = 0; lvl < numLevels; lvl++) {
      const EBLevelGrid&       eblg = m_petscGrid->getMFLevelGrids()[lvl]->getEBLevelGrid(iphase);
      const DisjointBoxLayout& dbl  = eblg.getDBL();

      const LevelData<BaseFab<PetscAMRCell>>& amrToPetsc = *m_petscGrid->getAMRToPetsc()[lvl];

      m_derivStencils[iphase][lvl] = RefCountedPtr<LayoutData<BaseIVFAB<PetscStencil>>>(
        new LayoutData<BaseIVFAB<PetscStencil>>(dbl));
      m_derivWeights[iphase][lvl] = RefCountedPtr<LayoutData<BaseIVFAB<Real>>>(new LayoutData<BaseIVFAB<Real>>(dbl));

      const DataIterator& dit  = dbl.dataIterator();
      const int           nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(+ : numWithoutStencil)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex&             din     = dit[mybox];
        const Box                    box     = dbl[din];
        const EBISBox&               ebisbox = eblg.getEBISL()[din];
        const EBGraph&               ebgraph = ebisbox.getEBGraph();
        const IntVectSet&            ivs     = ebisbox.getIrregIVS(box);
        const BaseFab<PetscAMRCell>& cells   = amrToPetsc[din];

        BaseIVFAB<PetscStencil>& stencils = (*m_derivStencils[iphase][lvl])[din];
        BaseIVFAB<Real>&         weights  = (*m_derivWeights[iphase][lvl])[din];

        stencils.define(ivs, ebgraph, 1);
        weights.define(ivs, ebgraph, 1);

        // Cut cells covered by a finer level have no row, and no stencil.
        for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit) {
          const VolIndex& vof = vofit();

          weights(vof, 0) = 0.0;
          stencils(vof, 0).clear();

          if (cells(vof.gridIndex()).getPetscRow(iphase) >= 0) {
            this->computeNormalDerivativeStencil(weights(vof, 0), stencils(vof, 0), iphase, lvl, din, vof);

            if (weights(vof, 0) == 0.0 && ebisbox.normal(vof) != RealVect::Zero) {
              numWithoutStencil++;
            }
          }
        }
      }
    }
  }

  m_numCellsWithoutStencil = ParallelOps::sum(numWithoutStencil);
  m_isDefined              = true;
}

void
PetscHelmholtzEBStencils::getNormalDerivativeStencil(Real&            a_boundaryWeight,
                                                     PetscStencil&    a_stencil,
                                                     const int        a_phase,
                                                     const int        a_level,
                                                     const DataIndex& a_din,
                                                     const VolIndex&  a_vof) const noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::getNormalDerivativeStencil");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBStencils::getNormalDerivativeStencil" << endl;
  }

  CH_assert(m_isDefined);
  CH_assert(a_phase >= 0 && a_phase < m_petscGrid->getNumPhases());
  CH_assert(a_level >= 0 && a_level < m_derivStencils[a_phase].size());

  a_boundaryWeight = (*m_derivWeights[a_phase][a_level])[a_din](a_vof, 0);
  a_stencil        = (*m_derivStencils[a_phase][a_level])[a_din](a_vof, 0);
}

void
PetscHelmholtzEBStencils::computeNormalDerivativeStencil(Real&            a_boundaryWeight,
                                                         PetscStencil&    a_stencil,
                                                         const int        a_phase,
                                                         const int        a_level,
                                                         const DataIndex& a_din,
                                                         const VolIndex&  a_vof) const noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::computeNormalDerivativeStencil");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBStencils::computeNormalDerivativeStencil" << endl;
  }

  const EBISBox& ebisbox = m_petscGrid->getMFLevelGrids()[a_level]->getEBLevelGrid(a_phase).getEBISL()[a_din];

  a_boundaryWeight = 0.0;
  a_stencil.clear();

  // Without a normal there is no EB face, and no flux through it.
  if (ebisbox.normal(a_vof) == RealVect::Zero) {
    return;
  }

  const VofUtils::Neighborhood neighborhoods[3] = {VofUtils::Neighborhood::SemiCircle,
                                                   VofUtils::Neighborhood::Quadrant,
                                                   VofUtils::Neighborhood::Radius};

  for (const VofUtils::Neighborhood neighborhood : neighborhoods) {
    for (int order = m_order; order > 0; order--) {
      if (this->getLeastSquaresStencil(a_boundaryWeight,
                                       a_stencil,
                                       a_phase,
                                       a_level,
                                       a_din,
                                       a_vof,
                                       neighborhood,
                                       order)) {
        return;
      }
    }
  }

  // No stencil was found. As in FieldSolverGMG, the cell then has no flux through its EB face: w_B = 0 and an empty
  // stencil.
  a_boundaryWeight = 0.0;
  a_stencil.clear();
}

bool
PetscHelmholtzEBStencils::getLeastSquaresStencil(Real&                        a_boundaryWeight,
                                                 PetscStencil&                a_stencil,
                                                 const int                    a_phase,
                                                 const int                    a_level,
                                                 const DataIndex&             a_din,
                                                 const VolIndex&              a_vof,
                                                 const VofUtils::Neighborhood a_neighborhood,
                                                 const int                    a_order) const noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::getLeastSquaresStencil");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBStencils::getLeastSquaresStencil" << endl;
  }

  constexpr bool addStartVof = false;
  constexpr Real deltaThresh = 0.5;

  const int  finestLevel = m_petscGrid->getMFLevelGrids().size() - 1;
  const bool hasFine     = a_level < finestLevel;
  const Real dx          = m_dx[a_level];
  const Real dxFine      = hasFine ? m_dx[a_level + 1] : 0.0;
  const int  refRat      = hasFine ? m_petscGrid->getRefinementRatios()[a_level] : 1;

  const EBISBox& ebisbox = m_petscGrid->getMFLevelGrids()[a_level]->getEBLevelGrid(a_phase).getEBISL()[a_din];
  const RealVect normal  = ebisbox.normal(a_vof);

  const BaseFab<PetscAMRCell>& cells = (*m_petscGrid->getAMRToPetsc()[a_level])[a_din];

  // The neighbourhood on this level, as in LeastSquares::getBndryGradSten.
  Vector<VolIndex> vofs;

  switch (a_neighborhood) {
  case VofUtils::Neighborhood::Quadrant: {
    vofs = VofUtils::getVofsInQuadrant(a_vof,
                                       ebisbox,
                                       normal,
                                       a_order,
                                       VofUtils::Connectivity::MonotonePath,
                                       addStartVof);

    break;
  }
  case VofUtils::Neighborhood::Radius: {
    vofs = VofUtils::getVofsInRadius(a_vof, ebisbox, a_order, VofUtils::Connectivity::MonotonePath, addStartVof);

    break;
  }
  case VofUtils::Neighborhood::SemiCircle: {
    vofs = VofUtils::getVofsInSemiCircle(a_vof,
                                         ebisbox,
                                         normal,
                                         a_order,
                                         deltaThresh,
                                         VofUtils::Connectivity::MonotonePath,
                                         Location::Cell::Boundary,
                                         Location::Cell::Center,
                                         addStartVof);

    break;
  }
  default: {
    MayDay::Abort("PetscHelmholtzEBStencils::getLeastSquaresStencil -- logic bust");

    break;
  }
  }

  // Sort the cells: rows and coarse-fine ghost cells stay on this level, and cells covered by the finer level are
  // replaced by the finer cells under them that are rows. Anything else, e.g. a cell outside the PETSc map, is not
  // used.
  Vector<VolIndex> vofsFine;
  Vector<VolIndex> vofsThis;

  for (int i = 0; i < vofs.size(); i++) {
    const IntVect& iv = vofs[i].gridIndex();

    if (!(cells.box().contains(iv))) {
      continue;
    }

    const PetscAMRCell& cell = cells(iv);

    if (cell.getPetscRow(a_phase) >= 0 || (cell.isGhostCF() && a_level > 0)) {
      vofsThis.push_back(vofs[i]);
    }
    else if (cell.isCoveredByFinerGrid() && hasFine) {
      // The refined grids of this level (FiCo) share this level's DataIndex, so a_din indexes them. Their map holds the
      // finer level's rows within its ghost cells; finer cells under a covered cell further out than that are not
      // used, like any cell outside the map on this level.
      const BaseFab<PetscAMRCell>& cellsFine = (*m_petscGrid->getAMRToPetscFiCo()[a_level + 1])[a_din];

      const Box fineBox = refine(Box(iv, iv), refRat);

      for (BoxIterator bit(fineBox); bit.ok(); ++bit) {
        const IntVect& ivFine = bit();

        if (cellsFine.box().contains(ivFine) && cellsFine(ivFine).getPetscRow(a_phase) >= 0) {
          vofsFine.push_back(VolIndex(ivFine, 0));
        }
      }
    }
  }

  // As in FieldSolverGMG, every cell of the neighbourhood enters the fit.
  const int numUnknowns = LeastSquares::getTaylorExpansionSize(a_order) - 1;

  if (static_cast<int>(vofsThis.size() + vofsFine.size()) <= numUnknowns) {
    return false;
  }

  // Displacements from the expansion point to the cell centers.
  const EBISBox&
    ebisboxFine = hasFine ? m_petscGrid->getMFLevelGridsFiCo()[a_level + 1]->getEBLevelGrid(a_phase).getEBISL()[a_din]
                          : ebisbox;

  Vector<RealVect> displacementsThis;
  Vector<RealVect> displacementsFine;

  const RealVect shift = this->getExpansionShift(a_phase, a_level, a_din, a_vof);

  for (int i = 0; i < vofsThis.size(); i++) {
    displacementsThis.push_back(
      LeastSquares::displacement(Location::Cell::Boundary, Location::Cell::Center, a_vof, vofsThis[i], ebisbox, dx) -
      shift);
  }

  for (int i = 0; i < vofsFine.size(); i++) {
    displacementsFine.push_back(LeastSquares::displacement(Location::Cell::Boundary,
                                                           Location::Cell::Center,
                                                           a_vof,
                                                           vofsFine[i],
                                                           ebisbox,
                                                           ebisboxFine,
                                                           dx,
                                                           dxFine) -
                                shift);
  }

  // Fit the gradient with phi_B known. The first stencil of each pair is over this level, the second over the finer.
  IntVectSet derivs;
  IntVectSet knownTerms;

  for (int dir = 0; dir < SpaceDim; dir++) {
    derivs |= BASISV(dir);
  }

  knownTerms |= IntVect::Zero;

  const auto gradient = LeastSquares::computeDualLevelStencils<Real>(derivs,
                                                                     knownTerms,
                                                                     vofsThis,
                                                                     vofsFine,
                                                                     displacementsThis,
                                                                     displacementsFine,
                                                                     m_weight,
                                                                     a_order);

  // Project onto n = -normal. The fit eliminated phi_B, which enters with minus the sum of the other weights.
  VoFStencil stencilThis;
  VoFStencil stencilFine;

  for (int dir = 0; dir < SpaceDim; dir++) {
    VoFStencil derivThis = gradient.at(BASISV(dir)).first;
    VoFStencil derivFine = gradient.at(BASISV(dir)).second;

    derivThis *= -normal[dir];
    derivFine *= -normal[dir];

    stencilThis += derivThis;
    stencilFine += derivFine;
  }

  // An empty result means that the least-squares solve failed.
  if (stencilThis.size() + stencilFine.size() == 0) {
    return false;
  }

  a_boundaryWeight = -(LeastSquares::sumAllWeights(stencilThis) + LeastSquares::sumAllWeights(stencilFine));

  // The normal derivative must increase with the boundary value. As in FieldSolverGMG, a fit with w_B <= 0, e.g. from
  // cells that lie nearly in the plane of the boundary, is rejected so that another neighbourhood or a lower order is
  // tried. On an interface, w_B <= 0 on one side can make b_0 w_0 + b_1 w_1 vanish.
  if (a_boundaryWeight <= 0.0) {
    a_boundaryWeight = 0.0;

    return false;
  }

  const RefCountedPtr<EBMultigridInterpolator> interpolator = (a_level > 0)
                                                                ? m_interpolators[a_level].getInterpolator(a_phase)
                                                                : RefCountedPtr<EBMultigridInterpolator>();

  a_stencil = m_petscGrid->getPetscStencilWithGhostCells(VoFStencil(),
                                                         stencilThis,
                                                         stencilFine,
                                                         interpolator,
                                                         a_level,
                                                         a_din,
                                                         a_phase);

  return true;
}

RealVect
PetscHelmholtzEBStencils::getBoundaryPosition(const int        a_phase,
                                              const int        a_level,
                                              const DataIndex& a_din,
                                              const VolIndex&  a_vof) const noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::getBoundaryPosition");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBStencils::getBoundaryPosition" << endl;
  }

  CH_assert(m_isDefined);

  const EBISBox& ebisbox = m_petscGrid->getMFLevelGrids()[a_level]->getEBLevelGrid(a_phase).getEBISL()[a_din];

  return m_probLo + Location::position(Location::Cell::Boundary, a_vof, ebisbox, m_dx[a_level]);
}

long long
PetscHelmholtzEBStencils::getNumCellsWithoutStencil() const noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::getNumCellsWithoutStencil");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBStencils::getNumCellsWithoutStencil" << endl;
  }

  CH_assert(m_isDefined);

  return m_numCellsWithoutStencil;
}

RealVect
PetscHelmholtzEBStencils::getExpansionShift(const int        a_phase,
                                            const int        a_level,
                                            const DataIndex& a_din,
                                            const VolIndex&  a_vof) const noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::getExpansionShift");

  RealVect shift = RealVect::Zero;

  if (!(m_interface.isNull())) {
    constexpr int  maxIter  = 10;
    constexpr Real maxShift = 0.5;

    const EBISBox& ebisbox = m_petscGrid->getMFLevelGrids()[a_level]->getEBLevelGrid(a_phase).getEBISL()[a_din];

    const Real     dx = m_dx[a_level];
    const RealVect xB = m_probLo + Location::position(Location::Cell::Boundary, a_vof, ebisbox, dx);
    const RealVect xS = PolyUtils::projectOntoSurface(m_interface, xB, 1.E-3 * dx, 1.E-10 * dx, maxIter);

    shift = xS - xB;

    // A point half a cell or more away is not on the part of the boundary that cuts this cell.
    if (!(shift.vectorLength() < maxShift * dx)) {
      shift = RealVect::Zero;
    }
  }

  return shift;
}

#include <CD_NamespaceFooter.H>

#endif
