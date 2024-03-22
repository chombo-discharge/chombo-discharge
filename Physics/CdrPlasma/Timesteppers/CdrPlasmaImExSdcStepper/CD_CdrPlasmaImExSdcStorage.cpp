/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaImExSdcStorage.cpp
  @brief  Implementation of CD_CdrPlasmaImExSdcStorage.H
  @author Robert Marskar
*/

// Our includes
#include <CD_CdrPlasmaImExSdcStepper.H>
#include <CD_CdrPlasmaImExSdcStorage.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaImExSdcStepper::CdrStorage::CdrStorage()
{}

CdrPlasmaImExSdcStepper::CdrStorage::CdrStorage(const RefCountedPtr<AmrMesh>& a_amr,
                                                const std::string             a_realm,
                                                const phase::which_phase      a_phase,
                                                const int                     a_ncomp)
{
  m_amr   = a_amr;
  m_realm = a_realm;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

CdrPlasmaImExSdcStepper::CdrStorage::~CdrStorage()
{
  deallocateStorage();
}

void
CdrPlasmaImExSdcStepper::CdrStorage::allocateStorage(const int a_p)
{
  m_p = a_p;

  m_amr->allocate(m_scratch, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratch2, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_error, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_old, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_divF, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_gradient, m_realm, m_phase, SpaceDim);
  m_amr->allocate(m_scratchD, m_realm, m_phase, SpaceDim);

  m_amr->allocate(m_scratchIVs, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratchIVD, m_realm, m_phase, SpaceDim);

  m_amr->allocate(m_scratchIV1, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV2, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV3, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV4, m_realm, m_phase, m_ncomp);

  m_amr->allocate(m_scratchIF1, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF2, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF3, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF4, m_realm, m_phase, m_ncomp);

  m_phi.resize(1 + m_p);
  m_FAR.resize(1 + m_p);
  m_FD.resize(1 + m_p);
  m_F.resize(1 + m_p);

  for (int m = 0; m <= m_p; m++) {
    m_amr->allocate(m_phi[m], m_realm, m_phase, m_ncomp);
    m_amr->allocate(m_FAR[m], m_realm, m_phase, m_ncomp);
    m_amr->allocate(m_FD[m], m_realm, m_phase, m_ncomp);
    m_amr->allocate(m_F[m], m_realm, m_phase, m_ncomp);
  }
}

void
CdrPlasmaImExSdcStepper::CdrStorage::deallocateStorage()
{

  m_amr->deallocate(m_scratch);
  m_amr->deallocate(m_scratch2);
  m_amr->deallocate(m_error);
  m_amr->deallocate(m_old);
  m_amr->deallocate(m_divF);
  m_amr->deallocate(m_gradient);
  m_amr->deallocate(m_scratchD);

  m_amr->deallocate(m_scratchIVs);
  m_amr->deallocate(m_scratchIVD);
  m_amr->deallocate(m_scratchIV1);
  m_amr->deallocate(m_scratchIV2);
  m_amr->deallocate(m_scratchIV3);
  m_amr->deallocate(m_scratchIV4);

  m_amr->deallocate(m_scratchIF1);
  m_amr->deallocate(m_scratchIF2);
  m_amr->deallocate(m_scratchIF3);
  m_amr->deallocate(m_scratchIF4);

  for (int m = 0; m <= m_p; m++) {
    m_amr->deallocate(m_phi[m]);
    m_amr->deallocate(m_FAR[m]);
    m_amr->deallocate(m_FD[m]);
    m_amr->deallocate(m_F[m]);
  }
}

CdrPlasmaImExSdcStepper::FieldStorage::FieldStorage()
{}

CdrPlasmaImExSdcStepper::FieldStorage::FieldStorage(const RefCountedPtr<AmrMesh>& a_amr,
                                                    const std::string             a_realm,
                                                    const phase::which_phase      a_phase,
                                                    const int                     a_ncomp)
{
  m_amr   = a_amr;
  m_realm = a_realm;
  m_ncomp = a_ncomp;
  m_phase = a_phase;
}

CdrPlasmaImExSdcStepper::FieldStorage::~FieldStorage()
{
  deallocateStorage();
}

void
CdrPlasmaImExSdcStepper::FieldStorage::allocateStorage(const int a_p)
{
  m_p = a_p;

  m_amr->allocate(m_previous, m_realm, m_ncomp);
  m_amr->allocate(m_E_cell, m_realm, m_phase, SpaceDim);
  m_amr->allocate(m_E_face, m_realm, m_phase, SpaceDim);
  m_amr->allocate(m_E_eb, m_realm, m_phase, SpaceDim);
  m_amr->allocate(m_E_dom, m_realm, m_phase, SpaceDim);
}

void
CdrPlasmaImExSdcStepper::FieldStorage::deallocateStorage()
{
  m_amr->deallocate(m_previous);
  m_amr->deallocate(m_E_cell);
  m_amr->deallocate(m_E_face);
  m_amr->deallocate(m_E_eb);
  m_amr->deallocate(m_E_dom);
}

CdrPlasmaImExSdcStepper::RtStorage::RtStorage()
{}

CdrPlasmaImExSdcStepper::RtStorage::RtStorage(const RefCountedPtr<AmrMesh>& a_amr,
                                              const std::string             a_realm,
                                              const phase::which_phase      a_phase,
                                              const int                     a_ncomp)
{
  m_amr   = a_amr;
  m_realm = a_realm;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

CdrPlasmaImExSdcStepper::RtStorage::~RtStorage()
{
  deallocateStorage();
}

void
CdrPlasmaImExSdcStepper::RtStorage::allocateStorage(const int a_p)
{
  m_p = a_p;

  m_amr->allocate(m_previous, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF, m_realm, m_phase, m_ncomp);
}

void
CdrPlasmaImExSdcStepper::RtStorage::deallocateStorage()
{
  m_amr->deallocate(m_previous);
  m_amr->deallocate(m_scratchIV);
  m_amr->deallocate(m_scratchIF);
}

CdrPlasmaImExSdcStepper::SigmaStorage::SigmaStorage()
{}

CdrPlasmaImExSdcStepper::SigmaStorage::SigmaStorage(const RefCountedPtr<AmrMesh>& a_amr,
                                                    const std::string             a_realm,
                                                    const phase::which_phase      a_phase,
                                                    const int                     a_ncomp)
{
  m_amr   = a_amr;
  m_realm = a_realm;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

CdrPlasmaImExSdcStepper::SigmaStorage::~SigmaStorage()
{
  deallocateStorage();
}

void
CdrPlasmaImExSdcStepper::SigmaStorage::allocateStorage(const int a_p)
{
  m_p = a_p;

  m_amr->allocate(m_scratch, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_error, m_realm, m_phase, m_ncomp);

  m_sigma.resize(1 + m_p);
  m_Fnew.resize(1 + m_p);
  m_Fold.resize(1 + m_p);

  for (int m = 0; m <= m_p; m++) {
    m_amr->allocate(m_sigma[m], m_realm, m_phase, m_ncomp);
    m_amr->allocate(m_Fnew[m], m_realm, m_phase, m_ncomp);
    m_amr->allocate(m_Fold[m], m_realm, m_phase, m_ncomp);
  }
}

void
CdrPlasmaImExSdcStepper::SigmaStorage::deallocateStorage()
{
  m_amr->deallocate(m_scratch);
  m_amr->deallocate(m_error);

  for (int m = 0; m <= m_p; m++) {
    m_amr->deallocate(m_sigma[m]);
    m_amr->deallocate(m_Fnew[m]);
    m_amr->deallocate(m_Fold[m]);
  }
}

#include <CD_NamespaceFooter.H>
