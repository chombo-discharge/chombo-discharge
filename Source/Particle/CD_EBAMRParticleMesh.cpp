/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBAMRParticleMesh.cpp
  @brief  Implementation of CD_EBAMRParticleMesh.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBAMRParticleMesh.H>
#include <CD_NamespaceHeader.H>

EBAMRParticleMesh::EBAMRParticleMesh()
{
  CH_TIME("EBAMRParticleMesh::EBAMRParticleMesh()");

  m_isDefined = false;
}

EBAMRParticleMesh::EBAMRParticleMesh(const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgs,
                                     const Vector<int>&                        a_refRat,
                                     const Vector<Real>&                       a_dx,
                                     const RealVect&                           a_probLo,
                                     const IntVect&                            a_ghost,
                                     const int                                 a_maxParticleWidth,
                                     const int                                 a_finestLevel)
{
  CH_TIME("EBAMRParticleMesh::EBAMRParticleMesh(full)");

  this->define(a_eblgs, a_refRat, a_dx, a_probLo, a_ghost, a_maxParticleWidth, a_finestLevel);
}

EBAMRParticleMesh::~EBAMRParticleMesh() { CH_TIME("EBAMRParticleMesh::~EBAMRParticleMesh()"); }

void
EBAMRParticleMesh::define(const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgs,
                          const Vector<int>&                        a_refRat,
                          const Vector<Real>&                       a_dx,
                          const RealVect&                           a_probLo,
                          const IntVect&                            a_ghost,
                          const int                                 a_maxParticleWidth,
                          const int                                 a_finestLevel)
{
  CH_TIME("EBAMRParticleMesh::define");

  m_eblgs            = a_eblgs;
  m_refRat           = a_refRat;
  m_dx               = a_dx;
  m_probLo           = a_probLo;
  m_ghost            = a_ghost;
  m_maxParticleWidth = a_maxParticleWidth;
  m_finestLevel      = a_finestLevel;

  this->defineLevelMotion();
  this->defineCoarseFineMotion();
  this->defineEBParticleMesh();

  m_isDefined = true;
}

void
EBAMRParticleMesh::defineLevelMotion()
{
  CH_TIME("EBAMRParticleMesh::defineLevelMotion");

  // TLDR: Define level Copiers. These are defined such that we can move data from valid+ghost -> valid. We need this because when we deposit particles
  //       we will also deposit into ghost cells that overlap with patches on the same level. The data in those ghost cells needs to find its way into
  //       the neighboring patches, which is what these Copiers do.
  m_levelCopiers.resize(1 + m_finestLevel);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const EBLevelGrid&       eblg   = *m_eblgs[lvl];
    const ProblemDomain&     domain = eblg.getDomain();
    const DisjointBoxLayout& dbl    = eblg.getDBL();

    // Define Copier. Note that the below code is basically the same as ghostDefine().
    const bool doExchange = true;

    // Define Copier as going from valid -> valid+ghost.
    m_levelCopiers[lvl].define(dbl, dbl, domain, m_ghost, doExchange);

    // Define Copier as going from valid+ghost -> valid.
    m_levelCopiers[lvl].reverse();
  }
}

void
EBAMRParticleMesh::defineCoarseFineMotion()
{
  CH_TIME("EBAMRParticleMesh::defineCoarseFineMotion");

  m_coarseFinePM.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {

    const bool hasCoar = (lvl > 0);

    if (hasCoar) {
      m_coarseFinePM[lvl] = RefCountedPtr<EBCoarseFineParticleMesh>(
        new EBCoarseFineParticleMesh(*m_eblgs[lvl - 1], *m_eblgs[lvl], m_refRat[lvl - 1], m_ghost));
    }
    else {
      m_coarseFinePM[lvl] = RefCountedPtr<EBCoarseFineParticleMesh>(nullptr);
    }
  }
}

void
EBAMRParticleMesh::defineEBParticleMesh()
{
  CH_TIME("EBAMRParticleMesh::defineEBParticleMesh");

  m_ebParticleMesh.resize(1 + m_finestLevel);
  m_ebParticleMeshFiCo.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const ProblemDomain&     domain = m_eblgs[lvl]->getDomain();
    const DisjointBoxLayout& dbl    = m_eblgs[lvl]->getDBL();
    const EBISLayout&        ebisl  = m_eblgs[lvl]->getEBISL();

    const bool hasCoar = lvl > 0;

    m_ebParticleMesh[lvl] = RefCountedPtr<LayoutData<EBParticleMesh>>(new LayoutData<EBParticleMesh>(dbl));

    // Define the "regular" particle-mesh interpolation objects. These are defined on the input grids.
    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box      cellBox = dbl[dit()];
      const EBISBox& ebisBox = ebisl[dit()];

      EBParticleMesh& particleMesh = (*m_ebParticleMesh[lvl])[dit()];

      particleMesh.define(domain, cellBox, ebisBox, m_dx[lvl] * RealVect::Unit, m_probLo);
    }

    // These are "special" particle-mesh interpolation objects for when we need to deposit coarse-level particles on a refined grid.
    if (hasCoar) {
      const EBLevelGrid& eblgFiCo = m_coarseFinePM[lvl]->getEblgFiCo();

      const ProblemDomain&     domain    = eblgFiCo.getDomain();
      const DisjointBoxLayout& dblFiCo   = eblgFiCo.getDBL();
      const EBISLayout&        ebislFiCo = eblgFiCo.getEBISL();

      m_ebParticleMeshFiCo[lvl] = RefCountedPtr<LayoutData<EBParticleMesh>>(new LayoutData<EBParticleMesh>(dblFiCo));

      for (DataIterator dit(dblFiCo); dit.ok(); ++dit) {
        const Box      cellBox = dblFiCo[dit()];
        const EBISBox& ebisBox = ebislFiCo[dit()];

        EBParticleMesh& particleMesh = (*m_ebParticleMeshFiCo[lvl])[dit()];

        particleMesh.define(domain, cellBox, ebisBox, m_dx[lvl] * RealVect::Unit, m_probLo);
      }
    }
    else {
      m_ebParticleMeshFiCo[lvl] = RefCountedPtr<LayoutData<EBParticleMesh>>(nullptr);
    }
  }
}

const EBParticleMesh&
EBAMRParticleMesh::getEBParticleMesh(const int a_lvl, const DataIndex& a_dit) const
{
  CH_TIME("EBAMRParticleMesh::getEBParticleMesh");

  CH_assert(a_lvl >= 0);
  CH_assert(a_lvl <= m_finestLevel);

  return (*m_ebParticleMesh[a_lvl])[a_dit];
}

#include <CD_NamespaceFooter.H>
