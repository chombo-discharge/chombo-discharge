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
                                   const Vector<RefCountedPtr<EBMultigridInterpolator>>&    a_coarseFineInterpolators,
                                   const Vector<int>&                                       a_refinementRatios,
                                   const Vector<Real>&                                      a_dx,
                                   const int&                                               a_finestLevel,
                                   const int&                                               a_ebbcOrder,
                                   const int&                                               a_ebbcWeight,
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
               a_ebbcOrder,
               a_ebbcWeight,
               a_jumpOrder,
               a_jumpWeight,
               a_dataLocation);
}

MFHelmholtzPetsc::~MFHelmholtzPetsc() noexcept
{
  CH_TIME("MFHelmholtzPetsc::~MFHelmholtzPetsc");
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
                         const Vector<RefCountedPtr<EBMultigridInterpolator>>&    a_coarseFineInterpolators,
                         const Vector<int>&                                       a_refinementRatios,
                         const Vector<Real>&                                      a_dx,
                         const int&                                               a_finestLevel,
                         const int&                                               a_ebbcOrder,
                         const int&                                               a_ebbcWeight,
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
  m_interpolators       = a_coarseFineInterpolators;
  m_refinementRatios    = a_refinementRatios;
  m_dx                  = a_dx;
  m_finestLevel         = a_finestLevel;
  m_ebbcOrder           = a_ebbcOrder;
  m_ebbcWeight          = a_ebbcWeight;
  m_jumpOrder           = a_jumpOrder;
  m_jumpWeight          = a_jumpWeight;
  m_dataLocation        = a_dataLocation;

  CH_assert(a_petscGrid->isDefined());
  CH_assert(m_finestLevel >= 0);
  CH_assert(m_ebbcOrder > 0);
  CH_assert(m_ebbcWeight >= 0);
  CH_assert(m_jumpOrder > 0);
  CH_assert(m_jumpWeight >= 0);

  m_isDefined = true;
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
MFHelmholtzPetsc::setEBBCOrder(const int a_order) noexcept
{
  CH_TIME("MFHelmholtzPetsc::setEBBCOrder");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::setEBBCOrder" << endl;
  }

  m_ebbcOrder = a_order;

  CH_assert(m_ebbcOrder > 0);
}

void
MFHelmholtzPetsc::setEBBCWeight(const int a_weight) noexcept
{
  CH_TIME("MFHelmholtzPetsc::setEBBCWeight");
  if (m_verbose) {
    pout() << "MFHelmholtzPetsc::setEBBCWeight" << endl;
  }

  m_ebbcWeight = a_weight;

  CH_assert(m_ebbcWeight >= 0);
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

#include <CD_NamespaceFooter.H>

#endif
