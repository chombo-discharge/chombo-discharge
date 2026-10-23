/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_RefinementBoxTagger.cpp
 * @brief  Implementation of CD_RefinementBoxTagger.H
 * @author Robert Marskar
 */

// Std includes
#include <algorithm>

// Chombo includes
#include <CH_Timer.H>
#include <BoxIterator.H>

// Our includes
#include <CD_RefinementBoxTagger.H>
#include <CD_ParallelOps.H>
#include <CD_NamespaceHeader.H>

RefinementBoxTagger::RefinementBoxTagger(const RefCountedPtr<AmrMesh>& a_amr, const std::string& a_realm) noexcept
  : m_amr(a_amr), m_realm(a_realm)
{
  CH_TIME("RefinementBoxTagger::RefinementBoxTagger");

  m_name = "RefinementBoxTagger";
}

RefinementBoxTagger::~RefinementBoxTagger() noexcept
{
  CH_TIME("RefinementBoxTagger::~RefinementBoxTagger");
}

void
RefinementBoxTagger::parseOptions()
{
  CH_TIME("RefinementBoxTagger::parseOptions");

  this->parseVerbosity();
  this->parseRefinementBoxes();
  this->parseBuffer();
}

void
RefinementBoxTagger::regrid()
{
  CH_TIME("RefinementBoxTagger::regrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::regrid()" << endl;
  }
}

bool
RefinementBoxTagger::tagCells(EBAMRTags& a_tags)
{
  CH_TIME("RefinementBoxTagger::tagCells");
  if (m_verbosity > 5) {
    pout() << m_name + "::tagCells(EBAMRTags)" << endl;
  }

  CH_assert(!m_amr.isNull());

  const RealVect probLo = m_amr->getProbLo();

  int haveTags = 0;

  // Tags on a level ask for the next finer level, so there are none on the finest allowed level.
  const int numTagLevels = std::min(static_cast<int>(a_tags.size()), m_amr->getMaxAmrDepth());

  for (int lvl = 0; lvl < numTagLevels; lvl++) {
    const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[lvl];
    const Real               dx   = m_amr->getDx()[lvl];
    const DataIterator&      dit  = dbl.dataIterator();
    const int                nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(max : haveTags)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din  = dit[mybox];
      DenseIntVectSet& tags = (*a_tags[lvl])[din];

      for (BoxIterator bit(dbl[din]); bit.ok(); ++bit) {
        const RealVect pos = probLo + (RealVect(bit()) + 0.5 * RealVect::Unit) * dx;

        if (this->getManualRefinementLevel(pos) > lvl) {
          tags |= bit();

          haveTags = 1;
        }
      }
    }
  }

  return ParallelOps::max(haveTags) > 0;
}

#include <CD_NamespaceFooter.H>
