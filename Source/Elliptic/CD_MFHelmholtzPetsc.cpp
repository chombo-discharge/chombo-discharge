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
#include <CD_VofUtils.H>
#include <CD_Timer.H>
#include <CD_MFHelmholtzPetsc.H>
#include <CD_NamespaceHeader.H>

#warning \
  "Before moving on with this class, fix the indexing logic in PetscGrid. Each rank must be able to reach into global rows, not just local ones"
#warning "Work item #1: Compute face flux stencils (interior)"

MFHelmholtzPetsc::MFHelmholtzPetsc() noexcept
{
  CH_TIME("MFHelmholtzPetsc::MFHelmholtzPetsc");

  m_isDefined = false;
  m_profile   = false;
  m_debug     = false;
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
{
  CH_TIME("MFHelmholtzPetsc::MFHelmholtzPetsc(full)");

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

std::pair<PetscScalar, PetscStencil>
MFHelmholtzPetsc::computeDirichletEBGradStencil(const VolIndex&  a_vof,
                                                const int&       a_phase,
                                                const int&       a_level,
                                                const DataIndex& a_din) const noexcept
{
  CH_TIME("MFHelmholtzPetsc::computeDirichletEBGradStencil");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::computeDirichletEBGradStencil" << endl;
  }

#warning "I need to expose the various coarse/fine stuff, too"
  const Vector<RefCountedPtr<MFLevelGrid>>& grids     = m_petscGrid->getMFLevelGrids();
  const Vector<RefCountedPtr<MFLevelGrid>>& gridsFiCo = m_petscGrid->getMFLevelGridsFiCo();
  const Vector<RefCountedPtr<MFLevelGrid>>& gridsCoFi = m_petscGrid->getMFLevelGridsCoFi();
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
