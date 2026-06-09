/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_RegularGeometry.cpp
  @brief  Implementation of CD_RegularGeometry.H
  @author Robert Marskar
*/

// Our includes
#include <CD_RegularGeometry.H>
#include <CD_NamespaceHeader.H>

RegularGeometry::RegularGeometry()
{
  m_dielectrics.resize(0);
  m_electrodes.resize(0);
}

RegularGeometry::~RegularGeometry()
{}

#include <CD_NamespaceFooter.H>
