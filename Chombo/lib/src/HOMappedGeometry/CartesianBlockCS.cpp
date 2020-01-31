#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CartesianBlockCS.H"
#include "BoxIterator.H"
#include "FourthOrderUtil.H"

#include "NamespaceHeader.H"

CartesianBlockCS::CartesianBlockCS(int a_blockNum,
                                   const IntVect& a_blockIndex,
                                   const RealVect& a_dx,
                                   const Box& a_blockBox)
{
  m_blockNum = a_blockNum;
  m_blockIndex = a_blockIndex;
  m_dx = a_dx;
  m_blockBox = a_blockBox;
}

CartesianBlockCS::~CartesianBlockCS()
{
}

RealVect
CartesianBlockCS::realCoord(const RealVect& a_Xi) const
{
  RealVect x = a_Xi - RealVect(m_blockIndex);
  return x;
}

RealVect
CartesianBlockCS::mappedCoord(const RealVect& a_x) const
{
  RealVect xi = a_x + RealVect(m_blockIndex);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      switch (m_blockIndex[idir])
        {
        case 0: // valid x[idir] in 0:1, valid a_Xi[idir] in 0:1
          if (a_x[idir] > 1.5)
            { // wrap around
              xi[idir] -= 2.;
            }
          break;
        case 1: // valid x[idir] in 1:2, valid a_Xi[idir] in 2:3
          if (a_x[idir] < 0.5)
            { // wrap around
              xi[idir] += 2.;
            }
          break;
        }
    }
  return xi;
}

Real
CartesianBlockCS::dXdXi(const RealVect& a_Xi,
                        int a_dirX,
                        int a_dirXi) const
{
  Real deriv = (a_dirX == a_dirXi) ? 1. : 0.;
  return deriv;
}

void
CartesianBlockCS::volFlux(FluxBox& a_volFlux,
                          const FluxBox& a_Nt,
                          const Box& a_box) const
{
  // This is adapted from CubedSphere2DPanelCS::volFlux().

  Box bx1 = grow(a_box, 1);
  CH_assert(a_Nt.box().contains(bx1));
  FluxBox NtX(bx1, SpaceDim);

  // FIXME: this assumes it is part of DoubleCartesianCS?
  // When m_blockIndex[idir] == 0, signs[idir] == 1;
  // when m_blockIndex[idir] == 1, signs[idir] == -1.
  IntVect signs = IntVect::Unit - 2 * m_blockIndex;

  // Midpoints are at 0.5 and 2.5.
  RealVect midXi = RealVect(2*m_blockIndex) + 0.5*RealVect::Unit;

  // Compute the point values of NtX.
  NtX.setVal(0.0);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      RealVect offset = 0.5*m_dx;
      offset[dir] = 0.0;

      FArrayBox& thisNtXdir = NtX[dir];
      BoxIterator bit(thisNtXdir.box());
      // this is going to be slow, but we can
      // eventually move this into fortran
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect mappedLoc = m_dx*iv;
          mappedLoc += offset;
          for (int comp=m_volInterval.begin(); comp<=m_volInterval.end(); comp++)
            {
              // only use components in volInterval
              thisNtXdir(iv, comp) = signs[comp] *
                (mappedLoc[comp] - midXi[comp]) / Real(m_volInterval.size());
            }
        }
    } // end loop over directions

  // Convert point values to 4th-order face averages
  fourthOrderAverageFace(NtX);
  // a_volFlux[dir] = <X[dir]>
  for (int dir=0; dir<SpaceDim; dir++)
    {
      // a_volFlux[dir][0] := NtX[dir][dir]
      a_volFlux[dir].copy(NtX[dir], dir, 0);
    }

  // Compute the volume flux on the faces.
  // computeMetricTermProductAverage(a_volFlux, X, a_Nt, SpaceDim, X, a_box);
}

#include "NamespaceFooter.H"
