/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_EBParticleMesh.cpp
  @brief  Implementationof CD_EBParticleMesh.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <CH_Timer.H>

// Our includes
#include <CD_EBParticleMesh.H>
#include <CD_NamespaceHeader.H>

EBParticleMesh::EBParticleMesh()
{
  CH_TIME("EBParticleMesh::EBParticleMesh");

  m_verbose = false;
}

EBParticleMesh::EBParticleMesh(const ProblemDomain& a_domain,
                               const Box&           a_region,
                               const EBISBox&       a_ebisbox,
                               const RealVect&      a_dx,
                               const RealVect&      a_probLo)
{
  CH_TIME("EBParticleMesh::EBParticleMesh");

  this->define(a_domain, a_region, a_ebisbox, a_dx, a_probLo);
}

void
EBParticleMesh::define(const ProblemDomain& a_domain,
                       const Box&           a_region,
                       const EBISBox&       a_ebisbox,
                       const RealVect&      a_dx,
                       const RealVect&      a_probLo)
{
  CH_TIME("EBParticleMesh::define");

  m_verbose = false;
  m_domain  = a_domain;
  m_region  = a_region;
  m_ebisbox = a_ebisbox;
  m_probLo  = a_probLo;
  m_dx      = a_dx;

  CH_assert(m_domain.contains(m_region));

  ParmParse pp("EBParticleMesh");
  pp.query("verbose", m_verbose);
}

#include <CD_NamespaceFooter.H>
