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

MFHelmholtzPetsc::MFHelmholtzPetsc(const RefCountedPtr<PetscGrid>&                       a_petscGrid,
                                   const RealVect&                                       a_probLo,
                                   const Real&                                           a_alpha,
                                   const Real&                                           a_beta,
                                   const Vector<RefCountedPtr<LevelData<MFCellFAB>>>&    a_aCoef,
                                   const Vector<RefCountedPtr<LevelData<MFFluxFAB>>>&    a_bCoef,
                                   const Vector<RefCountedPtr<LevelData<MFBaseIVFAB>>>&  a_bCoefIrreg,
                                   const Vector<RefCountedPtr<EBMultigridInterpolator>>& a_coarseFineInterpolators,
                                   const Vector<int>&                                    a_refinementRatios,
                                   const Vector<Real>&                                   a_dx,
                                   const int&                                            a_finestLevel,
                                   const Location::Cell&                                 a_dataLocation) noexcept
{
  CH_TIME("MFHelmholtzPetsc::MFHelmholtzPetsc(full)");

  this->define(a_petscGrid,
               a_probLo,
               a_alpha,
               a_beta,
               a_aCoef,
               a_bCoef,
               a_bCoefIrreg,
               a_coarseFineInterpolators,
               a_refinementRatios,
               a_dx,
               a_finestLevel,
               a_dataLocation);
}

MFHelmholtzPetsc::~MFHelmholtzPetsc() noexcept
{
  CH_TIME("MFHelmholtzPetsc::~MFHelmholtzPetsc");
}

void
MFHelmholtzPetsc::define(const RefCountedPtr<PetscGrid>&                       a_petscGrid,
                         const RealVect&                                       a_probLo,
                         const Real&                                           a_alpha,
                         const Real&                                           a_beta,
                         const Vector<RefCountedPtr<LevelData<MFCellFAB>>>&    a_aCoef,
                         const Vector<RefCountedPtr<LevelData<MFFluxFAB>>>&    a_bCoef,
                         const Vector<RefCountedPtr<LevelData<MFBaseIVFAB>>>&  a_bCoefIrreg,
                         const Vector<RefCountedPtr<EBMultigridInterpolator>>& a_coarseFineInterpolators,
                         const Vector<int>&                                    a_refinementRatios,
                         const Vector<Real>&                                   a_dx,
                         const int&                                            a_finestLevel,
                         const Location::Cell&                                 a_dataLocation) noexcept
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

  m_petscGrid        = a_petscGrid;
  m_probLo           = a_probLo;
  m_alpha            = a_alpha;
  m_beta             = a_beta;
  m_aCoef            = a_aCoef;
  m_bCoef            = a_bCoef;
  m_bCoefIrreg       = a_bCoefIrreg;
  m_interpolators    = a_coarseFineInterpolators;
  m_refinementRatios = a_refinementRatios;
  m_dx               = a_dx;
  m_finestLevel      = a_finestLevel;
  m_dataLocation     = a_dataLocation;

  CH_assert(a_petscGrid->isDefined());
  m_isDefined = true;
}

#include <CD_NamespaceFooter.H>

#endif
