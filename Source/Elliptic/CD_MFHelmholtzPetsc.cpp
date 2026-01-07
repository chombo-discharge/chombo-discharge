/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzPetsc.cpp
  @brief  Implementationof CD_MFHelmholtzPetsc.H
  @author Robert Marskar
*/

#ifdef CH_USE_PETSC

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_LeastSquares.H>
#include <CD_Timer.H>
#include <CD_MFHelmholtzPetsc.H>
#include <CD_NamespaceHeader.H>

#warning "Reconsider the usage of separate order for Dirichlet and jump stencils (don't want to compute twice)"
#warning "Work item #1: Compute Dirichlet EB stencils"
#warning "Work item #2: Compute face flux stencils (interior)"
#warning "Work item #3: Compute EB gradient stencils"
#warning "Work item #4: Compute domain flux stencils?"
#warning "Rename computeDiricihletDomain to computeDomainDirichlet"

MFHelmholtzPetsc::MFHelmholtzPetsc() noexcept
{
  CH_TIME("MFHelmholtzPetsc::MFHelmholtzPetsc");

  m_isDefined = false;
  m_profile   = false;
  m_debug     = false;

  ParmParse pp("MFHelmholtzPetsc");

  pp.query("debug", m_debug);
  pp.query("profile", m_profile);
  pp.query("verbose", m_verbose);
}

MFHelmholtzPetsc::MFHelmholtzPetsc(const RefCountedPtr<PetscGrid>&                          a_petscGrid,
                                   const RealVect&                                          a_probLo,
                                   const Real&                                              a_alpha,
                                   const Real&                                              a_beta,
                                   const bool&                                              a_multRhoByVolFrac,
                                   const bool&                                              a_multSigmaByAreaFrac,
                                   const Vector<RefCountedPtr<LevelData<MFCellFAB>>>&       a_aCoef,
                                   const Vector<RefCountedPtr<LevelData<MFFluxFAB>>>&       a_bCoef,
                                   const Vector<RefCountedPtr<LevelData<MFBaseIVFAB>>>&     a_bCoefIrreg,
                                   const Vector<RefCountedPtr<LevelData<MFCellFAB>>>&       a_rho,
                                   const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real>>>>& a_sigma,
                                   const Vector<MFMultigridInterpolator>&                   a_coarseFineInterpolators,
                                   const Vector<int>&                                       a_refinementRatios,
                                   const Vector<Real>&                                      a_dx,
                                   const int&                                               a_finestLevel,
                                   const int&                                               a_dirichletEBBCOrder,
                                   const int&                                               a_dirichletEBBCWeight,
                                   const int&                                               a_robinEBBCOrder,
                                   const int&                                               a_robinEBBCWeight,
                                   const int&                                               a_jumpOrder,
                                   const int&                                               a_jumpWeight,
                                   const Location::Cell&                                    a_dataLocation) noexcept
  : MFHelmholtzPetsc()
{
  CH_TIME("MFHelmholtzPetsc::MFHelmholtzPetsc(full)");

  m_ebDirichletNeighborhoods = {VofUtils::Neighborhood::SemiCircle,
                                VofUtils::Neighborhood::Quadrant,
                                VofUtils::Neighborhood::Radius};
  m_ebRobinNeighborhoods     = {VofUtils::Neighborhood::SemiCircle,
                                VofUtils::Neighborhood::Quadrant,
                                VofUtils::Neighborhood::Radius};
  m_ebJumpNeighborhoods      = {VofUtils::Neighborhood::SemiCircle,
                                VofUtils::Neighborhood::Quadrant,
                                VofUtils::Neighborhood::Radius};

  this->define(a_petscGrid,
               a_probLo,
               a_alpha,
               a_beta,
               a_multRhoByVolFrac,
               a_multSigmaByAreaFrac,
               a_aCoef,
               a_bCoef,
               a_bCoefIrreg,
               a_rho,
               a_sigma,
               a_coarseFineInterpolators,
               a_refinementRatios,
               a_dx,
               a_finestLevel,
               a_dirichletEBBCOrder,
               a_dirichletEBBCWeight,
               a_robinEBBCOrder,
               a_robinEBBCWeight,
               a_jumpOrder,
               a_jumpWeight,
               a_dataLocation);
}

MFHelmholtzPetsc::~MFHelmholtzPetsc() noexcept
{
  CH_TIME("MFHelmholtzPetsc::~MFHelmholtzPetsc");
}

bool
MFHelmholtzPetsc::isDefined() const noexcept
{
  CH_TIME("MFHelmholtzPetsc::isDefined");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::isDefined" << endl;
  }

  return m_isDefined;
}

void
MFHelmholtzPetsc::define(const RefCountedPtr<PetscGrid>&                          a_petscGrid,
                         const RealVect&                                          a_probLo,
                         const Real&                                              a_alpha,
                         const Real&                                              a_beta,
                         const bool&                                              a_multRhoByVolFrac,
                         const bool&                                              a_multSigmaByAreaFrac,
                         const Vector<RefCountedPtr<LevelData<MFCellFAB>>>&       a_aCoef,
                         const Vector<RefCountedPtr<LevelData<MFFluxFAB>>>&       a_bCoef,
                         const Vector<RefCountedPtr<LevelData<MFBaseIVFAB>>>&     a_bCoefIrreg,
                         const Vector<RefCountedPtr<LevelData<MFCellFAB>>>&       a_rho,
                         const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real>>>>& a_sigma,
                         const Vector<MFMultigridInterpolator>&                   a_coarseFineInterpolators,
                         const Vector<int>&                                       a_refinementRatios,
                         const Vector<Real>&                                      a_dx,
                         const int&                                               a_finestLevel,
                         const int&                                               a_dirichletEBBCOrder,
                         const int&                                               a_dirichletEBBCWeight,
                         const int&                                               a_robinEBBCOrder,
                         const int&                                               a_robinEBBCWeight,
                         const int&                                               a_jumpOrder,
                         const int&                                               a_jumpWeight,
                         const Location::Cell&                                    a_dataLocation) noexcept
{
  CH_TIME("MFHelmholtzPetsc::define");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::define" << endl;
  }

  ParmParse pp("MFHelmholtzPetsc");

  m_isDefined = false;
  m_debug     = false;
  m_profile   = false;
  m_verbose   = false;

  pp.query("debug", m_debug);
  pp.query("profile", m_profile);
  pp.query("verbose", m_verbose);

  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::define" << endl;
  }

  m_petscGrid           = a_petscGrid;
  m_probLo              = a_probLo;
  m_alpha               = a_alpha;
  m_beta                = a_beta;
  m_multRhoByVolFrac    = a_multRhoByVolFrac;
  m_multSigmaByAreaFrac = a_multSigmaByAreaFrac;
  m_aCoef               = a_aCoef;
  m_bCoef               = a_bCoef;
  m_bCoefIrreg          = a_bCoefIrreg;
  m_rho                 = a_rho;
  m_sigma               = a_sigma;
  m_interpolators       = a_coarseFineInterpolators;
  m_refinementRatios    = a_refinementRatios;
  m_dx                  = a_dx;
  m_finestLevel         = a_finestLevel;
  m_dirichletEBBCOrder  = a_dirichletEBBCOrder;
  m_dirichletEBBCWeight = a_dirichletEBBCWeight;
  m_robinEBBCOrder      = a_robinEBBCOrder;
  m_robinEBBCWeight     = a_robinEBBCWeight;
  m_jumpOrder           = a_jumpOrder;
  m_jumpWeight          = a_jumpWeight;
  m_dataLocation        = a_dataLocation;

  const int numLevels = 1 + m_finestLevel;

  CH_assert(a_petscGrid->isDefined());
  CH_assert(m_finestLevel >= 0);
  CH_assert(m_dirichletEBBCOrder > 0);
  CH_assert(m_dirichletEBBCWeight >= 0);
  CH_assert(m_robinEBBCOrder > 0);
  CH_assert(m_robinEBBCWeight >= 0);
  CH_assert(m_jumpOrder > 0);
  CH_assert(m_jumpWeight >= 0);
  CH_assert(m_aCoef.size() == numLevels);
  CH_assert(m_bCoef.size() == numLevels);
  CH_assert(m_bCoefIrreg.size() == numLevels);
  CH_assert(m_rho.size() == numLevels);
  CH_assert(m_sigma.size() == numLevels);
  CH_assert(m_interpolators.size() == numLevels);
  CH_assert(m_dx.size() == numLevels);
  CH_assert(m_refinementRatios.size() == numLevels);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    CH_assert(m_aCoef[lvl]->isDefined());
    CH_assert(m_bCoef[lvl]->isDefined());
    CH_assert(m_bCoefIrreg[lvl]->isDefined());
    CH_assert(m_rho[lvl]->isDefined());
    CH_assert(m_sigma[lvl]->isDefined());
    CH_assert(m_refinementRatios[lvl] >= 2);
    CH_assert(m_refinementRatios[lvl] % 2 == 0);
  }

  m_numPhases = m_petscGrid->getNumPhases();
  m_isDefined = true;

#warning "debug code enabled -- users should be responsible for this one"
#if 1
  this->computeEBDirichletStencils();
  this->computeEBGradStencils();
#endif
}

void
MFHelmholtzPetsc::clear() noexcept
{
  CH_TIME("MFHelmholtzPetsc::clear");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::clear" << endl;
  }

  m_petscGrid = RefCountedPtr<PetscGrid>(nullptr);
  m_aCoef.resize(0);
  m_bCoef.resize(0);
  m_bCoefIrreg.resize(0);
  m_isDefined = false;
}

void
MFHelmholtzPetsc::setMultRhoByVolFrac(const bool a_multRhoByVolFrac) noexcept
{
  CH_TIME("MFHelmholtzPetsc::setMultRhoByVolFrac");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::setMultRhoByVolFrac" << endl;
  }

  m_multRhoByVolFrac = a_multRhoByVolFrac;
}

void
MFHelmholtzPetsc::setMultSigmaByAreaFrac(const bool a_multSigmaByAreaFrac) noexcept
{
  CH_TIME("MFHelmholtzPetsc::setMultSigmaByAreaFrac");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::setMultSigmaByAreaFrac" << endl;
  }

  m_multSigmaByAreaFrac = a_multSigmaByAreaFrac;
}

void
MFHelmholtzPetsc::setDirichletEBBCOrder(const int a_order) noexcept
{
  CH_TIME("MFHelmholtzPetsc::setDirichletEBBCOrder");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::setDirichletEBBCOrder" << endl;
  }

  m_dirichletEBBCOrder = a_order;

  CH_assert(m_dirichletEBBCOrder > 0);
}

void
MFHelmholtzPetsc::setDirichletEBBCWeight(const int a_weight) noexcept
{
  CH_TIME("MFHelmholtzPetsc::setDirichletEBBCWeight");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::setDirichletEBBCWeight" << endl;
  }

  m_dirichletEBBCWeight = a_weight;

  CH_assert(m_dirichletEBBCWeight >= 0);
}

void
MFHelmholtzPetsc::setRobinEBBCOrder(const int a_order) noexcept
{
  CH_TIME("MFHelmholtzPetsc::setRobinEBBCOrder");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::setRobinEBBCOrder" << endl;
  }

  m_robinEBBCOrder = a_order;

  CH_assert(m_robinEBBCOrder > 0);
}

void
MFHelmholtzPetsc::setRobinEBBCWeight(const int a_weight) noexcept
{
  CH_TIME("MFHelmholtzPetsc::setRobinEBBCWeight");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::setRobinEBBCWeight" << endl;
  }

  m_robinEBBCWeight = a_weight;

  CH_assert(m_robinEBBCWeight >= 0);
}

void
MFHelmholtzPetsc::setJumpOrder(const int a_order) noexcept
{
  CH_TIME("MFHelmholtzPetsc::setJumpOrder");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::setJumpOrder" << endl;
  }

  m_jumpOrder = a_order;

  CH_assert(m_jumpOrder > 0);
}

void
MFHelmholtzPetsc::setJumpWeight(const int a_weight) noexcept
{
  CH_TIME("MFHelmholtzPetsc::setJumpWeight");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::setJumpWeight" << endl;
  }

  m_jumpWeight = a_weight;

  CH_assert(m_jumpWeight >= 0);
}

void
MFHelmholtzPetsc::setEBDirichletNeighborhoods(const std::vector<VofUtils::Neighborhood>& a_neighborhoods) noexcept
{
  CH_TIME("MFHelmholtzPetsc::setEBDirichletNeighborhoods");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::setEBDirichletNeighborhoods" << endl;
  }

  CH_assert(a_neighborhoods.size() > 0);

  m_ebDirichletNeighborhoods = a_neighborhoods;
}

void
MFHelmholtzPetsc::setEBRobinNeighborhoods(const std::vector<VofUtils::Neighborhood>& a_neighborhoods) noexcept
{
  CH_TIME("MFHelmholtzPetsc::setEBRobinNeighborhoods");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::setEBRobinNeighborhoods" << endl;
  }

  CH_assert(a_neighborhoods.size() > 0);

  m_ebRobinNeighborhoods = a_neighborhoods;
}

void
MFHelmholtzPetsc::setEBJumpNeighborhoods(const std::vector<VofUtils::Neighborhood>& a_neighborhoods) noexcept
{
  CH_TIME("MFHelmholtzPetsc::setEBJumpNeighborhoods");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::setEBJumpNeighborhoods" << endl;
  }

  CH_assert(a_neighborhoods.size() > 0);

  m_ebJumpNeighborhoods = a_neighborhoods;
}

void
MFHelmholtzPetsc::computeEBDirichletStencils() noexcept
{
  CH_TIME("MFHelmholtzPetsc::computeEBDirichletStencils");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::computeEBDirichletStencils" << endl;
  }

  CH_assert(m_isDefined);

  constexpr int comp    = 0;
  constexpr int numComp = 1;

  for (int iphase = 0; iphase < m_numPhases; iphase++) {
    m_ebDirichletStencils[iphase].resize(1 + m_finestLevel);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const MFLevelGrid&       mflg  = *(m_petscGrid->getMFLevelGrids()[lvl]);
      const EBLevelGrid&       eblg  = mflg.getEBLevelGrid(iphase);
      const EBISLayout&        ebisl = eblg.getEBISL();
      const DisjointBoxLayout& dbl   = eblg.getDBL();
      const DataIterator&      dit   = dbl.dataIterator();
      const int                nbox  = dit.size();

      m_ebDirichletStencils[iphase][lvl] = RefCountedPtr<LayoutData<BaseIVFAB<std::pair<PetscScalar, PetscStencil>>>>(
        new LayoutData<BaseIVFAB<std::pair<PetscScalar, PetscStencil>>>(dbl));

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex&             din        = dit[mybox];
        const Box&                   cellBox    = dbl[din];
        const EBISBox&               ebisBox    = ebisl[din];
        const EBGraph&               ebgraph    = ebisBox.getEBGraph();
        const BaseFab<PetscAMRCell>& petscToAMR = (*(m_petscGrid->getAMRToPetsc()[lvl]))[din];
        const IntVectSet&            irregIVS   = (*m_petscGrid->getIrregIVS(iphase)[lvl])[din];

        VoFIterator&                                     vofit    = (*m_petscGrid->getVoFIterator(iphase)[lvl])[din];
        BaseIVFAB<std::pair<PetscScalar, PetscStencil>>& stencils = (*m_ebDirichletStencils[iphase][lvl])[din];

        stencils.define(irregIVS, ebgraph, numComp);

        auto kernel = [&](const VolIndex& vof) -> void {
          const PetscAMRCell petscCell = petscToAMR(vof.gridIndex(), comp);
          const PetscInt     petscRow  = petscCell.getPetscRow(iphase);

          if (petscRow < 0) {
            MayDay::Abort("MFHelmholtzPetsc::computeEBDirichletStencils -- logic bust 'petscRow < 0'");
          }

          stencils(vof, comp) = this->computeEBDirichletStencil(vof, iphase, lvl, din);
        };

        BoxLoops::loop(vofit, kernel);
      }
    }
  }
}

std::pair<PetscScalar, PetscStencil>
MFHelmholtzPetsc::computeEBDirichletStencil(const VolIndex&  a_vof,
                                            const int&       a_phase,
                                            const int&       a_level,
                                            const DataIndex& a_din) const noexcept
{
  CH_TIME("MFHelmholtzPetsc::computeEBDirichletStencil");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::computeEBDirichletStencil" << endl;
  }

  CH_assert(a_phase < m_numPhases);

  const bool hasFine = a_level < m_finestLevel;

  Real dx;
  Real dxFine;

  MFLevelGrid mflg;
  MFLevelGrid mflgFine;

  EBLevelGrid eblg;
  EBLevelGrid eblgFine;

  EBISLayout ebisl;
  EBISLayout ebislFine;

  EBISBox ebisBox;
  EBISBox ebisBoxFine;

  Vector<VolIndex> vofs;
  Vector<VolIndex> vofsFine;

  mflg    = *m_petscGrid->getMFLevelGrids()[a_level];
  eblg    = mflg.getEBLevelGrid(a_phase);
  ebisl   = eblg.getEBISL();
  ebisBox = ebisl[a_din];
  dx      = m_dx[a_level];

  vofs.resize(0);
  vofsFine.resize(0);

  if (hasFine) {
    mflgFine    = *m_petscGrid->getMFLevelGridsFiCo()[a_level + 1];
    eblgFine    = mflg.getEBLevelGrid(a_phase);
    ebislFine   = eblgFine.getEBISL();
    ebisBoxFine = ebislFine[a_din];
    dxFine      = m_dx[a_level + 1];

    // First, we fetch all the vofs on the current level as usual.
    const Vector<VolIndex> allVoFs = VofUtils::getVofsInRadius(a_vof,
                                                               ebisBox,
                                                               m_dirichletEBBCOrder,
                                                               VofUtils::Connectivity::MonotonePath,
                                                               false);

    // Iterate through the vofs on this level. If they are covered by a finer grid we discard them, and instead
    // use the vofs defined by refinement of the current vof.
    for (int i = 0; i < allVoFs.size(); i++) {
      const Vector<VolIndex>& refinedVoFs = ebisBox.refine(vofs[i]);
    }
  }
  else {
    vofs = VofUtils::getVofsInRadius(a_vof, ebisBox, m_dirichletEBBCOrder, VofUtils::Connectivity::MonotonePath, false);
  }

  return std::make_pair(0, PetscStencil());
}

void
MFHelmholtzPetsc::computeEBGradStencils() noexcept
{
  CH_TIME("MFHelmholtzPetsc::computeEBGradStencils");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::computeEBGradStencils" << endl;
  }

  CH_assert(m_isDefined);

#warning "Continue development here"

  for (int iphase = 0; iphase < m_numPhases; iphase++) {
    m_ebGradStencils[iphase].resize(1 + m_finestLevel);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const MFLevelGrid&       mflg  = *(m_petscGrid->getMFLevelGrids()[lvl]);
      const EBLevelGrid&       eblg  = mflg.getEBLevelGrid(iphase);
      const EBISLayout&        ebisl = eblg.getEBISL();
      const DisjointBoxLayout& dbl   = eblg.getDBL();
      const DataIterator&      dit   = dbl.dataIterator();
      const int                nbox  = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex&  din     = dit[mybox];
        const Box&        cellBox = dbl[din];
        const EBISBox&    ebisBox = ebisl[din];
        const EBGraph&    ebgraph = ebisBox.getEBGraph();
        const IntVectSet& ivs     = ebisBox.getIrregIVS(cellBox);

        VoFIterator vofit(ivs, ebgraph);

        auto kernel = [&](const VolIndex& vof) -> void {

        };

        BoxLoops::loop(vofit, kernel);
      }
    }
  }
};

void
MFHelmholtzPetsc::computeEBExtrapStencils() noexcept
{
  CH_TIME("MFHelmholtzPetsc::computeEBExtrapStencils");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::computeEBExtrapStencils" << endl;
  }
}

PetscStencil
MFHelmholtzPetsc::computeInteriorFaceFluxStencil(const VolIndex&      a_vof,
                                                 const int&           a_phase,
                                                 const int&           a_level,
                                                 const int            a_dir,
                                                 const Side::LoHiSide a_side,
                                                 const DataIndex&     a_din) const noexcept
{
  CH_TIME("MFHelmholtzPetsc::computeInteriorFaceFluxStencil(EB)");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::computeInteriorFaceFluxStencil(EB)" << endl;
  }
}

PetscStencil
MFHelmholtzPetsc::computeInteriorFaceFluxStencil(const IntVect&       a_cell,
                                                 const int            a_phase,
                                                 const int            a_level,
                                                 const int            a_dir,
                                                 const Side::LoHiSide a_side,
                                                 const DataIndex&     a_din) const noexcept
{
  CH_TIME("MFHelmholtzPetsc::computeInteriorFaceFluxStencil(Regular)");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::computeInteriorFaceFluxStencil(Regular)" << endl;
  }
}

PetscStencil
MFHelmholtzPetsc::computeRobinEBGradStencil(const VolIndex&  a_vof,
                                            const int&       a_phase,
                                            const int&       a_level,
                                            const DataIndex& a_din) const noexcept
{
  CH_TIME("MFHelmholtzPetsc::computeRobinEBGradStencil");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::computeRobinEBGradStencil" << endl;
  }
}

PetscScalar
MFHelmholtzPetsc::computeNeumannEBGradStencil(const VolIndex&  a_vof,
                                              const int&       a_phase,
                                              const int&       a_level,
                                              const DataIndex& a_din) const noexcept
{
  CH_TIME("MFHelmholtzPetsc::computeNeumannEBGradStencil");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::computeNeumannEBGradStencil" << endl;
  }
}

std::pair<MFHelmholtzPetsc::PetscColumn, MFHelmholtzPetsc::PetscColumn>
MFHelmholtzPetsc::computeJumpFluxStencils(const IntVect&   a_cell,
                                          const int&       a_phase,
                                          const int&       a_level,
                                          const DataIndex& a_din) const noexcept
{
  CH_TIME("MFHelmholtzPetsc::computeJumpFluxStencils");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::computeJumpFluxStencils" << endl;
  }
}

std::pair<PetscStencil, PetscScalar>
MFHelmholtzPetsc::computeDirichletDomainGradStencil(const VolIndex&      a_vof,
                                                    const int&           a_phase,
                                                    const int&           a_level,
                                                    const int&           a_dir,
                                                    const Side::LoHiSide a_side,
                                                    const DataIndex&     a_din) const noexcept
{
  CH_TIME("MFHelmholtzPetsc::computeDirichletDomainGradStencil");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::computeDirichletDomainGradStencil" << endl;
  }
}

std::pair<PetscStencil, PetscScalar>
MFHelmholtzPetsc::computeRobinDomainGradStencil(const VolIndex&      a_vof,
                                                const int&           a_phase,
                                                const int&           a_level,
                                                const int&           a_dir,
                                                const Side::LoHiSide a_side,
                                                const DataIndex&     a_din) const noexcept
{
  CH_TIME("MFHelmholtzPetsc::computeRobinDomainGradStencil");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::computeRobinDomainGradStencil" << endl;
  }
}

PetscScalar
MFHelmholtzPetsc::computeNeumannDomainGradStencil(const VolIndex&      a_vof,
                                                  const int&           a_phase,
                                                  const int&           a_level,
                                                  const int&           a_dir,
                                                  const Side::LoHiSide a_side,
                                                  const DataIndex&     a_din) const noexcept
{
  CH_TIME("MFHelmholtzPetsc::computeNeumannDomainGradStencil");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::computeNeumannDomainGradStencil" << endl;
  }
}

MFHelmholtzPetsc::PetscColumn
MFHelmholtzPetsc::computeStencil(const VolIndex&  a_vof,
                                 const int&       a_phase,
                                 const int&       a_level,
                                 const DataIndex& a_din) const noexcept
{
  CH_TIME("MFHelmholtzPetsc::computeStencil");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::computeStencil" << endl;
  }
}

#include <CD_NamespaceFooter.H>

#endif
