/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ParticleMergeTagger.cpp
 * @brief  Implementation of CD_ParticleMergeTagger.H
 * @author Robert Marskar
 */

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_ParticleMergeTagger.H>
#include <CD_NamespaceHeader.H>

ParticleMergeTagger::ParticleMergeTagger() noexcept
{
  m_name      = "ParticleMergeTagger";
  m_verbosity = -1;
  m_refineDir = SpaceDim - 1;
  m_refineLo  = 0.0;
  m_refineHi  = 0.0;
}

ParticleMergeTagger::~ParticleMergeTagger() noexcept = default;

void
ParticleMergeTagger::define(const RefCountedPtr<AmrMesh>& a_amr) noexcept
{
  CH_TIME("ParticleMergeTagger::define");

  m_amr = a_amr;
}

void
ParticleMergeTagger::parseOptions()
{
  CH_TIME("ParticleMergeTagger::parseOptions");

  ParmParse pp(m_name.c_str());

  pp.get("refine_dir", m_refineDir);
  pp.get("refine_lo", m_refineLo);
  pp.get("refine_hi", m_refineHi);

  if (m_refineDir < 0 || m_refineDir >= SpaceDim) {
    MayDay::Error("ParticleMergeTagger::parseOptions - 'refine_dir' must be a coordinate direction");
  }
  if (m_refineHi <= m_refineLo) {
    MayDay::Error("ParticleMergeTagger::parseOptions - 'refine_hi' must exceed 'refine_lo'");
  }
}

void
ParticleMergeTagger::regrid()
{
  CH_TIME("ParticleMergeTagger::regrid");
}

bool
ParticleMergeTagger::tagCells(EBAMRTags& a_tags)
{
  CH_TIME("ParticleMergeTagger::tagCells");

  const RealVect probLo = m_amr->getProbLo();

  // Tag on every level that is allowed to be refined. The slab is a fixed physical region, so the tags
  // are the same set of cells every time -- what changes between levels is only the resolution the slab
  // is expressed at.
  for (int lvl = 0; lvl < a_tags.size(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_amr->getRealms()[0])[lvl];
    const Real               dx  = m_amr->getDx()[lvl];

    const DataIterator& dit  = dbl.dataIterator();
    const int           nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      DenseIntVectSet& tags = (*a_tags[lvl])[din];

      tags.makeEmptyBits();

      for (BoxIterator bit(dbl[din]); bit.ok(); ++bit) {
        const IntVect iv = bit();
        const Real    x  = probLo[m_refineDir] + dx * (iv[m_refineDir] + 0.5);

        if (x >= m_refineLo && x <= m_refineHi) {
          tags |= iv;
        }
      }
    }
  }

  return true;
}

void
ParticleMergeTagger::gatherTags(Vector<IntVectSet>& a_levelTags, const EBAMRTags& a_tags) const noexcept
{
  CH_TIME("ParticleMergeTagger::gatherTags");

  a_levelTags.resize(a_tags.size(), IntVectSet());

  for (int lvl = 0; lvl < a_tags.size(); lvl++) {
    IntVectSet& levelTags = a_levelTags[lvl];

    levelTags.makeEmpty();

    const DataIterator& dit  = a_tags[lvl]->dataIterator();
    const int           nbox = dit.size();

    // Serial on purpose: IntVectSet is not thread-safe to accumulate into, and this runs once per regrid
    // over a tag set that is a slab, not a field.
    for (int mybox = 0; mybox < nbox; mybox++) {
      levelTags |= IntVectSet((*a_tags[lvl])[dit[mybox]]);
    }
  }
}

#include <CD_NamespaceFooter.H>
