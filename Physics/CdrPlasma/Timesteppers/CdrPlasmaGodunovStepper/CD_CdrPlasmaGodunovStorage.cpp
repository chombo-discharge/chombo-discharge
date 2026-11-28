/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_CdrPlasmaGodunovStorage.cpp
 * @brief  Implementation of CD_CdrPlasmaGodunovStorage.H
 * @author Robert Marskar
 */

// Our includes
#include <CD_CdrPlasmaGodunovStepper.H>
#include <CD_CdrPlasmaGodunovStorage.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaGodunovStepper::CdrStorage::CdrStorage(const RefCountedPtr<AmrMesh>& a_amr,
                                                const std::string&            a_realm,
                                                const phase::which_phase      a_phase)
  : m_scratch(a_amr->fetchScratch<EBAMRCellData>(a_realm, a_phase, 1, a_amr->getNumberOfGhostCells())),
    m_scratch2(a_amr->fetchScratch<EBAMRCellData>(a_realm, a_phase, 1, a_amr->getNumberOfGhostCells())),
    m_cellExtr(a_amr->fetchScratch<EBAMRCellData>(a_realm, a_phase, 1, a_amr->getNumberOfGhostCells())),
    m_gradient(a_amr->fetchScratch<EBAMRCellData>(a_realm, a_phase, SpaceDim, a_amr->getNumberOfGhostCells())),
    m_scratchIV1(a_amr->fetchScratch<EBAMRIVData>(a_realm, a_phase, 1, a_amr->getNumberOfGhostCells())),
    m_scratchIV2(a_amr->fetchScratch<EBAMRIVData>(a_realm, a_phase, 1, a_amr->getNumberOfGhostCells())),
    m_scratchIV3(a_amr->fetchScratch<EBAMRIVData>(a_realm, a_phase, 1, a_amr->getNumberOfGhostCells())),
    m_scratchIV4(a_amr->fetchScratch<EBAMRIVData>(a_realm, a_phase, 1, a_amr->getNumberOfGhostCells())),
    m_scratchIF1(a_amr->fetchScratch<EBAMRIFData>(a_realm, a_phase, 1, a_amr->getNumberOfGhostCells())),
    m_scratchIF2(a_amr->fetchScratch<EBAMRIFData>(a_realm, a_phase, 1, a_amr->getNumberOfGhostCells())),
    m_scratchIF3(a_amr->fetchScratch<EBAMRIFData>(a_realm, a_phase, 1, a_amr->getNumberOfGhostCells())),
    m_scratchIF4(a_amr->fetchScratch<EBAMRIFData>(a_realm, a_phase, 1, a_amr->getNumberOfGhostCells()))
{}

CdrPlasmaGodunovStepper::CdrStorage::~CdrStorage() = default;

CdrPlasmaGodunovStepper::FieldStorage::FieldStorage(const RefCountedPtr<AmrMesh>& a_amr,
                                                    const std::string&            a_realm,
                                                    const phase::which_phase      a_phase)
  : m_electricFieldCell(a_amr->fetchScratch<EBAMRCellData>(a_realm, a_phase, SpaceDim, a_amr->getNumberOfGhostCells())),
    m_electricFieldEB(a_amr->fetchScratch<EBAMRIVData>(a_realm, a_phase, SpaceDim, a_amr->getNumberOfGhostCells())),
    m_electricFieldDomain(a_amr->fetchScratch<EBAMRIFData>(a_realm, a_phase, SpaceDim, a_amr->getNumberOfGhostCells()))
{}

CdrPlasmaGodunovStepper::FieldStorage::~FieldStorage() = default;

CdrPlasmaGodunovStepper::RtStorage::RtStorage(const RefCountedPtr<AmrMesh>& a_amr,
                                              const std::string&            a_realm,
                                              const phase::which_phase      a_phase)
  : m_scratchIV(a_amr->fetchScratch<EBAMRIVData>(a_realm, a_phase, 1, a_amr->getNumberOfGhostCells())),
    m_scratchIF(a_amr->fetchScratch<EBAMRIFData>(a_realm, a_phase, 1, a_amr->getNumberOfGhostCells()))
{}

CdrPlasmaGodunovStepper::RtStorage::~RtStorage() = default;

#include <CD_NamespaceFooter.H>
