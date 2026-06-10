/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
   @file   CD_BrownianWalkerTagger.cpp
   @brief  Implementation of CD_BrownianWalkerTagger.H
   @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_BrownianWalkerTagger.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::BrownianWalker;

BrownianWalkerTagger::BrownianWalkerTagger(RefCountedPtr<ItoSolver>& a_solver, RefCountedPtr<AmrMesh>& a_amr)
  : m_solver(a_solver), m_amr(a_amr)
{
  CH_TIME("BrownianWalkerTagger::BrownianWalkerTagger");

  m_name      = "BrownianWalker";
  m_verbosity = -1;
}

BrownianWalkerTagger::~BrownianWalkerTagger()
{
  CH_TIME("BrownianWalkerTagger::~BrownianWalkerTagger");
}

void
BrownianWalkerTagger::regrid()
{
  CH_TIME("BrownianWalkerTagger::regrid");
}

void
BrownianWalkerTagger::parseOptions()
{
  CH_TIME("BrownianWalkerTagger::parseOptions");

  ParmParse pp(m_name.c_str());
  pp.get("refine_magn", m_refMagn);

  this->parseBuffer();
}

bool
BrownianWalkerTagger::tagCells(EBAMRTags& /*a_tags*/)
{
  return true;
}

#include <CD_NamespaceFooter.H>
