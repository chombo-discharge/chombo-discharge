/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CoarseInterpQuadCF.cpp
  @brief  Implementation of CD_CoarseInterpQuadCF.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <DenseIntVectSet.H>
#include <NeighborIterator.H>

// Our includes
#include <CD_CoarseInterpQuadCF.H>
#include <CD_NamespaceHeader.H>

CoarseInterpQuadCF::CoarseInterpQuadCF() noexcept
{
  CH_TIME("CoarseInterpQuadCF::CoarseInterpQuadCF()");

  m_isDefined = false;
}

CoarseInterpQuadCF::~CoarseInterpQuadCF() noexcept { CH_TIME("CoarseInterpQuadCF::~CoarseInterpQuadCF()"); }

void
CoarseInterpQuadCF::define(const DisjointBoxLayout& a_dblFine,
                           const ProblemDomain&     a_domainCoar,
                           const DataIndex&         a_dit,
                           const Box&               a_fineGhostCells,
                           const int                a_refRat,
                           const int                a_ignoreDir) noexcept
{
  CH_TIME("CoarseInterpQuadCF::define");

  m_dblFine    = a_dblFine;
  m_domainCoar = a_domainCoar;
  m_dit        = a_dit;
  m_refRat     = a_refRat;
  m_ignoreDir  = a_ignoreDir;
  m_stencilBox = coarsen(a_fineGhostCells, m_refRat);

  // Define stencil objects.
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (dir != a_ignoreDir) {
      m_firstDerivStencils[dir].define(m_stencilBox, 1);
      m_secondDerivStencils[dir].define(m_stencilBox, 1);
    }
  }
#if CH_SPACEDIM == 3
  m_mixedDerivStencils.define(m_stencilBox, 1);
#endif

  // Compute stencils.
  this->defineFirstDerivStencils();
  this->defineSecondDerivStencils();
  this->defineMixedDerivStencils();

  m_isDefined = true;
}

void
CoarseInterpQuadCF::defineFirstDerivStencils() noexcept
{
  CH_TIME("CoarseInterpQuadCF::defineFirstDerivStencils");

  // Stencils will have a radius of 2, so figure out which cells we are allowed to use.
  DenseIntVectSet validCells(grow(m_stencilBox, 2), true);
  validCells -= m_dblFine[m_dit];
  NeighborIterator nit(m_dblFine);
  for (nit.begin(m_dit); nit.ok(); ++nit) {
    validCells -= m_dblFine[nit()];
  }
  validCells &= m_domainCoar;

  // Go through the cells and figure out which FD approximation we will use.
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (dir != m_ignoreDir) {
      for (BoxIterator bit(m_stencilBox); bit.ok(); ++bit) {
        const IntVect ivCoar = bit();

        const IntVect ivLoLo = ivCoar - 2 * BASISV(dir);
        const IntVect ivLo   = ivCoar - BASISV(dir);
        const IntVect ivHi   = ivCoar + BASISV(dir);
        const IntVect ivHiHi = ivCoar + 2 * BASISV(dir);

        const bool useLoLo = validCells[ivLoLo];
        const bool useLo   = validCells[ivLo];
        const bool useHi   = validCells[ivHi];
        const bool useHiHi = validCells[ivHiHi];

        // Switch between differencing formulas. We want to use second order differencing if we can, but use
        // first order if we must.
        if (useLo && useHi) {
          m_firstDerivStencils[dir](ivCoar, 0) = FirstDerivStencil::Centered2;
        }
        else if (useLo && useLoLo) {
          m_firstDerivStencils[dir](ivCoar, 0) = FirstDerivStencil::Backward2;
        }
        else if (useHi && useHiHi) {
          m_firstDerivStencils[dir](ivCoar, 0) = FirstDerivStencil::Forward2;
        }
        else if (useLo && !useHi) {
          m_firstDerivStencils[dir](ivCoar, 0) = FirstDerivStencil::Backward1;
        }
        else if (!useLo && useHi) {
          m_firstDerivStencils[dir](ivCoar, 0) = FirstDerivStencil::Forward1;
        }
      }
    }
  }
}

void
CoarseInterpQuadCF::defineSecondDerivStencils() noexcept
{
  CH_TIME("CoarseInterpQuadCF::defineSecondDerivStencils");

  // Stencils will have a radius of 2, so figure out which cells we are allowed to use.
  DenseIntVectSet validCells(grow(m_stencilBox, 2), true);
  validCells -= m_dblFine[m_dit];
  NeighborIterator nit(m_dblFine);
  for (nit.begin(m_dit); nit.ok(); ++nit) {
    validCells -= m_dblFine[nit()];
  }
  validCells &= m_domainCoar;

  // Go through the cells and compute finite difference approximations.
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (dir != m_ignoreDir) {
      BaseFab<DerivStencil>& derivStencilsDir = m_secondDerivStencils[dir];

      for (BoxIterator bit(m_stencilBox); bit.ok(); ++bit) {
        const IntVect ivCoar = bit();

        DerivStencil& derivSten = derivStencilsDir(ivCoar, 0);
        derivSten.clear();

        const IntVect ivLoLo = ivCoar - 2 * BASISV(dir);
        const IntVect ivLo   = ivCoar - BASISV(dir);
        const IntVect ivHi   = ivCoar + BASISV(dir);
        const IntVect ivHiHi = ivCoar + 2 * BASISV(dir);

        const bool useLoLo = validCells[ivLoLo];
        const bool useLo   = validCells[ivLo];
        const bool useHi   = validCells[ivHi];
        const bool useHiHi = validCells[ivHiHi];

        // Switch between differencing formulas. We want to use second order differencing if we can, but use
        // first order if we must.
        if (useLo && useHi) {
          // Centered differencing, which is second order accurate.
          derivSten.accumulate(ivLo, 1.0);
          derivSten.accumulate(ivCoar, -2.0);
          derivSten.accumulate(ivHi, 1.0);
        }
        else if (useLo && useLoLo) {
          // 2nd order backward differencing.
          derivSten.accumulate(ivLoLo, 1.0);
          derivSten.accumulate(ivLo, -2.0);
          derivSten.accumulate(ivCoar, 1.0);
        }
        else if (useHi && useHiHi) {
          // Forward differenceing
          derivSten.accumulate(ivCoar, 1.0);
          derivSten.accumulate(ivHi, -2.0);
          derivSten.accumulate(ivHiHi, 1.0);
        }
      }
    }
  }
}

void
CoarseInterpQuadCF::defineMixedDerivStencils() noexcept
{
#if CH_SPACEDIM == 3
  CH_TIME("CoarseInterpQuadCF::defineMixedDerivStencils");

  // Stencils will have a radius of 2, so figure out which cells we are allowed to use.
  DenseIntVectSet validCells(grow(m_stencilBox, 2), true);
  validCells -= m_dblFine[m_dit];
  NeighborIterator nit(m_dblFine);
  for (nit.begin(m_dit); nit.ok(); ++nit) {
    validCells -= m_dblFine[nit()];
  }
  validCells &= m_domainCoar;

  // Go through the cells and compute finite difference approximations to the mixed derivative. We do this ala Chombo
  // and average the mixed-derivative stencils on edges of the cell since some of the cells we otherwise would need
  // might be covered by the fine grid.
  for (BoxIterator bit(m_stencilBox); bit.ok(); ++bit) {
    const IntVect ivCoar = bit();

    DerivStencil& derivSten = m_mixedDerivStencils(ivCoar, 0);
    derivSten.clear();

    int dir1 = -1;
    int dir2 = -1;

    if (m_ignoreDir == 0) {
      dir1 = 1;
      dir2 = 2;
    }
    else if (m_ignoreDir == 1) {
      dir1 = 0;
      dir2 = 2;
    }
    else if (m_ignoreDir == 2) {
      dir1 = 0;
      dir2 = 1;
    }

    // Quadrant which may or may not be available.
    Vector<Box> quadrants;

    quadrants.push_back(Box(ivCoar, ivCoar + BASISV(dir1) + BASISV(dir2)));
    quadrants.push_back(Box(ivCoar - BASISV(dir1), ivCoar + BASISV(dir2)));
    quadrants.push_back(Box(ivCoar - BASISV(dir1) - BASISV(dir2), ivCoar));
    quadrants.push_back(Box(ivCoar - BASISV(dir2), ivCoar + BASISV(dir1)));

    int numQuadrants;

    for (int i = 0; i < quadrants.size(); i++) {
      const Box quad = quadrants[i];

      if (validCells.contains(quad)) {
        numQuadrants++;

        for (BoxIterator bit(quad); bit.ok(); ++bit) {
          derivSten.accumulate(bit(), 0.25);
        }
      }
    }

    if (numQuadrants > 0) {
      derivSten *= 1.0 / numQuadrants;
    }
  }
#endif
}

Real
CoarseInterpQuadCF::computeFirstDeriv(const FArrayBox& a_coarPhi,
                                      const IntVect&   a_ivCoar,
                                      const int        a_dir,
                                      const int        a_coarVar) const noexcept
{
  CH_TIME("CoarseInterpQuadCF::computeFirstDeriv");

  CH_assert(a_dir != m_ignoreDir);
  CH_assert(m_stencilBox.contains(a_ivCoar));

  Real firstDeriv = 0.0;

  const IntVect unitDir = BASISV(a_dir);

  switch (m_firstDerivStencils[a_dir](a_ivCoar, 0)) {
  case FirstDerivStencil::Centered2: {
    firstDeriv += 0.5 * a_coarPhi(a_ivCoar + unitDir, a_coarVar);
    firstDeriv -= 0.5 * a_coarPhi(a_ivCoar - unitDir, a_coarVar);

    break;
  }
  case FirstDerivStencil::Backward2: {
    firstDeriv += 1.5 * a_coarPhi(a_ivCoar, a_coarVar);
    firstDeriv -= 2.0 * a_coarPhi(a_ivCoar - unitDir, a_coarVar);
    firstDeriv += 0.5 * a_coarPhi(a_ivCoar - 2 * unitDir, a_coarVar);

    break;
  }
  case FirstDerivStencil::Forward2: {
    firstDeriv -= 1.5 * a_coarPhi(a_ivCoar, a_coarVar);
    firstDeriv += 2.0 * a_coarPhi(a_ivCoar + unitDir, a_coarVar);
    firstDeriv -= 0.5 * a_coarPhi(a_ivCoar + 2 * unitDir, a_coarVar);

    break;
  }
  case FirstDerivStencil::Backward1: {
    firstDeriv += a_coarPhi(a_ivCoar, a_coarVar);
    firstDeriv -= a_coarPhi(a_ivCoar - unitDir, a_coarVar);

    break;
  }
  case FirstDerivStencil::Forward1: {
    firstDeriv += a_coarPhi(a_ivCoar + unitDir, a_coarVar);
    firstDeriv -= a_coarPhi(a_ivCoar, a_coarVar);

    break;
  }
  default: {
    firstDeriv = 0.0;

    break;
  }
  }

  return firstDeriv;
}

Real
CoarseInterpQuadCF::computeSecondDeriv(const FArrayBox& a_coarPhi,
                                       const IntVect&   a_ivCoar,
                                       const int        a_dir,
                                       const int        a_coarVar) const noexcept
{
  CH_TIME("CoarseInterpQuadCF::computeSecondDeriv");

  CH_assert(a_dir != m_ignoreDir);
  CH_assert(m_stencilBox.contains(a_ivCoar));

  const DerivStencil& stencil = m_secondDerivStencils[a_dir](a_ivCoar, 0);

  Real secondDeriv = 0.0;
  for (int i = 0; i < stencil.size(); i++) {
    secondDeriv += stencil.getWeight(i) * a_coarPhi(stencil.getIndex(i), a_coarVar);
  }

  return secondDeriv;
}

Real
CoarseInterpQuadCF::computeMixedDeriv(const FArrayBox& a_coarPhi,
                                      const IntVect&   a_ivCoar,
                                      const int        a_coarVar) const noexcept
{
  CH_TIME("CoarseInterpQuadCF::computeMixedDeriv");

  CH_assert(m_stencilBox.contains(a_ivCoar));

  Real mixedDeriv = 0.0;
#if CH_SPACEDIM == 3
  const DerivStencil& stencil = m_mixedDerivStencils(a_ivCoar, 0);

  for (int i = 0; i < stencil.size(); i++) {
    mixedDeriv += stencil.getWeight(i) * a_coarPhi(stencil.getIndex(i), a_coarVar);
  }
#endif

  return mixedDeriv;
}

#include <CD_NamespaceFooter.H>
