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

  // Store the coordinates tangential to ignoreDir. Needed for mix-deriv stencil.
  if (m_ignoreDir == 0) {
    m_tanDir1 = 1;
    m_tanDir2 = 2;
  }
  else if (m_ignoreDir == 1) {
    m_tanDir1 = 0;
    m_tanDir2 = 2;
  }
  else if (m_ignoreDir == 2) {
    m_tanDir1 = 0;
    m_tanDir2 = 1;
  }

  // Define stencil objects.
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (dir != a_ignoreDir) {
      m_firstDerivStencils[dir].define(m_stencilBox, 1);
      m_secondDerivStencils[dir].define(m_stencilBox, 1);
    }
  }
#if CH_SPACEDIM == 3
  m_mixedDerivStencils.define(m_stencilBox, 1);
  m_explicitMixedDerivStencils.define(m_stencilBox, 1);
#endif

  // Compute stencils.
  this->defineStencils();
  this->defineMixedDerivStencils();

  m_isDefined = true;
}

void
CoarseInterpQuadCF::defineStencils() noexcept
{
  CH_TIME("CoarseInterpQuadCF::defineStencils");

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
          m_firstDerivStencils[dir](ivCoar, 0)  = FirstDerivStencil::Centered2;
          m_secondDerivStencils[dir](ivCoar, 0) = SecondDerivStencil::Centered2;
        }
        else if (useLo && useLoLo) {
          m_firstDerivStencils[dir](ivCoar, 0)  = FirstDerivStencil::Backward2;
          m_secondDerivStencils[dir](ivCoar, 0) = SecondDerivStencil::Backward1;
        }
        else if (useHi && useHiHi) {
          m_firstDerivStencils[dir](ivCoar, 0)  = FirstDerivStencil::Forward2;
          m_secondDerivStencils[dir](ivCoar, 0) = SecondDerivStencil::Forward1;
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

  constexpr Real factor = 1. / 16.;

  // Go through the cells and compute finite difference approximations to the mixed derivative. We do this ala Chombo
  // and average the mixed-derivative stencils on edges of the cell since some of the cells we otherwise would need
  // might be covered by the fine grid.
  for (BoxIterator bit(m_stencilBox); bit.ok(); ++bit) {
    const IntVect ivCoar = bit();

    if (validCells.contains(grow(Box(ivCoar, ivCoar), 1))) {
      m_mixedDerivStencils(ivCoar, 0) = MixedDerivStencil::Standard;
    }
    else {
      m_mixedDerivStencils(ivCoar, 0) = MixedDerivStencil::Explicit;

      DerivStencil& derivSten = m_explicitMixedDerivStencils(ivCoar, 0);
      derivSten.clear();

      // Quadrant which may or may not be available.
      int numQuadrants = 0;

      const Box baseQuadrant(ivCoar, ivCoar + BASISV(m_tanDir1) + BASISV(m_tanDir2));

      Vector<IntVect> shifts;
      shifts.push_back(IntVect::Zero);
      shifts.push_back(-BASISV(m_tanDir1));
      shifts.push_back(-BASISV(m_tanDir1) - BASISV(m_tanDir2));
      shifts.push_back(-BASISV(m_tanDir2));
      for (int i = 0; i < shifts.size(); i++) {
        Box shiftedBox = baseQuadrant;
        shiftedBox.shift(shifts[i]);

        if (validCells.contains(shiftedBox)) {
          numQuadrants++;

          // Standard mixed-derivative stencil. The true stencil is
          // d^2f/dxdx = (f(i+1,j+1) + f(i-1,j-1) - f(i+1,j-1) - f(i-1,j+1))/(4*dx*dx).
          // But we are computing on the edge between cells, so we're actually using dx = dxCoar/2, which is
          // why factor is 1./16 instead of 1/4
          derivSten.accumulate(ivCoar + shifts[i], factor);
          derivSten.accumulate(ivCoar + shifts[i] + BASISV(m_tanDir1) + BASISV(m_tanDir2), factor);
          derivSten.accumulate(ivCoar + shifts[i] + BASISV(m_tanDir1), -factor);
          derivSten.accumulate(ivCoar + shifts[i] + BASISV(m_tanDir2), -factor);
        }
      }

      // Do the average over all the mixed-deriv stencils we found in the available quadrants
      if (numQuadrants > 0) {
        derivSten *= 1.0 / numQuadrants;
      }
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

  Real secondDeriv = 0.0;

  const IntVect unitDir = BASISV(a_dir);

  switch (m_secondDerivStencils[a_dir](a_ivCoar, 0)) {
  case SecondDerivStencil::Centered2: {
    secondDeriv += a_coarPhi(a_ivCoar - unitDir, a_coarVar);
    secondDeriv -= 2. * a_coarPhi(a_ivCoar, a_coarVar);
    secondDeriv += a_coarPhi(a_ivCoar + unitDir, a_coarVar);

    break;
  }
  case SecondDerivStencil::Backward1: {
    secondDeriv += a_coarPhi(a_ivCoar - 2 * unitDir, a_coarVar);
    secondDeriv -= 2. * a_coarPhi(a_ivCoar - unitDir, a_coarVar);
    secondDeriv += a_coarPhi(a_ivCoar, a_coarVar);

    break;
  }
  case SecondDerivStencil::Forward1: {
    secondDeriv += a_coarPhi(a_ivCoar, a_coarVar);
    secondDeriv -= 2. * a_coarPhi(a_ivCoar + unitDir, a_coarVar);
    secondDeriv += a_coarPhi(a_ivCoar + 2 * unitDir, a_coarVar);

    break;
  }
  default: {
    secondDeriv = 0.0;

    break;
  }
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
  switch (m_mixedDerivStencils(a_ivCoar, 0)) {
  case MixedDerivStencil::Standard: {

    const IntVect tanDir1 = BASISV(m_tanDir1);
    const IntVect tanDir2 = BASISV(m_tanDir2);

    mixedDeriv += 0.25 * a_coarPhi(a_ivCoar + tanDir1 + tanDir2, a_coarVar);
    mixedDeriv += 0.25 * a_coarPhi(a_ivCoar - tanDir1 - tanDir2, a_coarVar);
    mixedDeriv += 0.25 * a_coarPhi(a_ivCoar + tanDir1 - tanDir2, a_coarVar);
    mixedDeriv += 0.25 * a_coarPhi(a_ivCoar - tanDir1 + tanDir2, a_coarVar);

    break;
  }
  case MixedDerivStencil::Explicit: {
    const DerivStencil& stencil = m_explicitMixedDerivStencils(a_ivCoar, 0);

    for (int i = 0; i < stencil.size(); i++) {
      mixedDeriv += stencil.getWeight(i) * a_coarPhi(stencil.getIndex(i), a_coarVar);
    }

    break;
  }
  default: {
    mixedDeriv = 0.0;

    break;
  }
  }

#endif

  return mixedDeriv;
}

#include <CD_NamespaceFooter.H>
