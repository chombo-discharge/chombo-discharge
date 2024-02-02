/* chombo-discharge
 * Copyright © 2024 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CellInfo.cpp
  @brief  Implementation of CD_CellInfo.H
  @author Robert Marskar
*/

// Our includes
#include <CD_DataOps.H>
#include <CD_CellInfo.H>
#include <CD_NamespaceHeader.H>

const CellInfo CellInfo::Regular = CellInfo();

CellInfo::CellInfo() noexcept
{
  m_volFrac       = 1.0;
  m_bndryCentroid = RealVect::Zero;
  m_bndryNormal   = RealVect::Zero;
  m_validLo       = -0.5 * RealVect::Unit;
  m_validHi       = 0.5 * RealVect::Unit;
}

CellInfo::CellInfo(const Real a_volFrac, const RealVect a_bndryCentroid, const RealVect a_bndryNormal) noexcept {
  m_volFrac = a_volFrac;
  m_bndryCentroid = a_bndryCentroid;
  m_bndryNormal = a_bndryNormal;
  
  DataOps::computeMinValidBox(m_validLo, m_validHi, m_bndryNormal, m_bndryCentroid);
}

CellInfo::~CellInfo() noexcept {}

Real&
CellInfo::getVolFrac() noexcept
{
  return (m_volFrac);
}

const Real&
CellInfo::getVolFrac() const noexcept
{
  return (m_volFrac);
}

RealVect&
CellInfo::getBndryCentroid() noexcept
{
  return (m_bndryCentroid);
}

const RealVect&
CellInfo::getBndryCentroid() const noexcept
{
  return (m_bndryCentroid);
}

RealVect&
CellInfo::getBndryNormal() noexcept
{
  return (m_bndryNormal);
}

const RealVect&
CellInfo::getBndryNormal() const noexcept
{
  return (m_bndryNormal);
}

RealVect&
CellInfo::getValidLo() noexcept
{
  return (m_validLo);
}

const RealVect&
CellInfo::getValidLo() const noexcept
{
  return (m_validLo);
}

RealVect&
CellInfo::getValidHi() noexcept
{
  return (m_validHi);
}

const RealVect&
CellInfo::getValidHi() const noexcept
{
  return (m_validHi);
}

#include <CD_NamespaceFooter.H>
