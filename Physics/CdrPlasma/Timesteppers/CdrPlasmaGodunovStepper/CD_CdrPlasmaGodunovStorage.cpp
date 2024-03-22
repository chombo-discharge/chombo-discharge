/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaGodunovStorage.cpp
  @brief  Implementation of CdrPlasmaGodunovStepper_storage.H
  @author Robert Marskar
  @date   Aug. 2019
*/

// Our includes
#include <CD_CdrPlasmaGodunovStepper.H>
#include <CD_CdrPlasmaGodunovStorage.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaGodunovStepper::CdrStorage::CdrStorage(const RefCountedPtr<AmrMesh>& a_amr,
                                                const std::string             a_realm,
                                                const phase::which_phase      a_phase)
{
  m_amr   = a_amr;
  m_realm = a_realm;
  m_phase = a_phase;
}

CdrPlasmaGodunovStepper::CdrStorage::~CdrStorage()
{
  deallocateStorage();
}

void
CdrPlasmaGodunovStepper::CdrStorage::allocateStorage()
{

  // Allocate all storage.
  constexpr int nComp = 1;

  // Allocate EBAMRCellData data holders
  m_amr->allocate(m_scratch, m_realm, m_phase, nComp);
  m_amr->allocate(m_scratch2, m_realm, m_phase, nComp);
  m_amr->allocate(m_cellExtr, m_realm, m_phase, nComp);
  m_amr->allocate(m_gradient, m_realm, m_phase, SpaceDim);

  // Allocate EBAMRIVData data holders
  m_amr->allocate(m_scratchIV1, m_realm, m_phase, nComp);
  m_amr->allocate(m_scratchIV2, m_realm, m_phase, nComp);
  m_amr->allocate(m_scratchIV3, m_realm, m_phase, nComp);
  m_amr->allocate(m_scratchIV4, m_realm, m_phase, nComp);

  // Allocate EBAMRIFData data holders
  m_amr->allocate(m_scratchIF1, m_realm, m_phase, nComp);
  m_amr->allocate(m_scratchIF2, m_realm, m_phase, nComp);
  m_amr->allocate(m_scratchIF3, m_realm, m_phase, nComp);
  m_amr->allocate(m_scratchIF4, m_realm, m_phase, nComp);
}

void
CdrPlasmaGodunovStepper::CdrStorage::deallocateStorage()
{

  // Deallocate EBAMRCellData storage
  m_amr->deallocate(m_scratch);
  m_amr->deallocate(m_scratch2);
  m_amr->deallocate(m_cellExtr);
  m_amr->deallocate(m_gradient);

  // Deallocate EBAMRIVData storage
  m_amr->deallocate(m_scratchIV1);
  m_amr->deallocate(m_scratchIV2);
  m_amr->deallocate(m_scratchIV3);
  m_amr->deallocate(m_scratchIV4);

  // Deallocate EBAMRIFData storage
  m_amr->deallocate(m_scratchIF1);
  m_amr->deallocate(m_scratchIF2);
  m_amr->deallocate(m_scratchIF3);
  m_amr->deallocate(m_scratchIF4);
}

CdrPlasmaGodunovStepper::FieldStorage::FieldStorage(const RefCountedPtr<AmrMesh>& a_amr,
                                                    const std::string             a_realm,
                                                    const phase::which_phase      a_phase)
{
  m_amr   = a_amr;
  m_phase = a_phase;
  m_realm = a_realm;
}

CdrPlasmaGodunovStepper::FieldStorage::~FieldStorage()
{
  this->deallocateStorage();
}

void
CdrPlasmaGodunovStepper::FieldStorage::allocateStorage()
{

  // Allocate the various data
  m_amr->allocate(m_electricFieldCell, m_realm, m_phase, SpaceDim);
  m_amr->allocate(m_electricFieldEB, m_realm, m_phase, SpaceDim);
  m_amr->allocate(m_electricFieldDomain, m_realm, m_phase, SpaceDim);
}

void
CdrPlasmaGodunovStepper::FieldStorage::deallocateStorage()
{

  // Deallocate the various data.
  m_amr->deallocate(m_electricFieldCell);
  m_amr->deallocate(m_electricFieldEB);
  m_amr->deallocate(m_electricFieldDomain);
}

CdrPlasmaGodunovStepper::RtStorage::RtStorage(const RefCountedPtr<AmrMesh>& a_amr,
                                              const std::string             a_realm,
                                              const phase::which_phase      a_phase)
{
  m_amr   = a_amr;
  m_realm = a_realm;
  m_phase = a_phase;
}

CdrPlasmaGodunovStepper::RtStorage::~RtStorage()
{}

void
CdrPlasmaGodunovStepper::RtStorage::allocateStorage()
{

  constexpr int nComp = 1;

  // Allocate the storage
  m_amr->allocate(m_scratchIV, m_realm, m_phase, nComp);
  m_amr->allocate(m_scratchIF, m_realm, m_phase, nComp);
}

void
CdrPlasmaGodunovStepper::RtStorage::deallocateStorage()
{

  // Deallocate the storage
  m_amr->deallocate(m_scratchIV);
  m_amr->deallocate(m_scratchIF);
}

CdrPlasmaGodunovStepper::SigmaStorage::SigmaStorage(const RefCountedPtr<AmrMesh>& a_amr,
                                                    const std::string             a_realm,
                                                    const phase::which_phase      a_phase)
{
  m_amr   = a_amr;
  m_realm = a_realm;
  m_phase = a_phase;
}

CdrPlasmaGodunovStepper::SigmaStorage::~SigmaStorage()
{}

void
CdrPlasmaGodunovStepper::SigmaStorage::allocateStorage()
{
  constexpr int nComp = 1;

  m_amr->allocate(m_scratch, m_realm, m_phase, nComp);
}

void
CdrPlasmaGodunovStepper::SigmaStorage::deallocateStorage()
{
  m_amr->deallocate(m_scratch);
}

#include <CD_NamespaceFooter.H>
