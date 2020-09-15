#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CartesianBlockRotateCS.H"
#include "BoxIterator.H"
#include "FourthOrderUtil.H"

#include "NamespaceHeader.H"

CartesianBlockRotateCS::CartesianBlockRotateCS(int a_blockNum,
                                               const IntVect& a_blockIndex,
                                               const RealVect& a_dx,
                                               const Box& a_blockBox,
                                               IntVect a_perm,
                                               IntVect a_sign)
{
  m_blockNum = a_blockNum;
  m_blockIndex = a_blockIndex;
  m_dx = a_dx;
  m_blockBox = a_blockBox;
  m_perm = a_perm;
  m_sign = a_sign;
  // Get m_blockCenter, with components either 0.5 or 2.5.
  IntVect shiftInt = 2 * m_blockIndex;
  m_blockCenter = RealVect(shiftInt) + 0.5 * RealVect::Unit;
  for (int dirSrc = 0; dirSrc < SpaceDim; dirSrc++)
    {
      IntVect ivSrc = BASISV(dirSrc);
      for (int dirDst = 0; dirDst < SpaceDim; dirDst++)
        {
          int ival = m_sign[dirDst] * ivSrc[m_perm[dirDst]]; // ivDst[dirDst]
          Real val = ival * 1.;
          m_vecTransform.setComponent(dirDst, dirSrc, val);
        }
    }
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}

CartesianBlockRotateCS::~CartesianBlockRotateCS()
{
}

RealVect
CartesianBlockRotateCS::realCoord(const RealVect& a_Xi) const
{
  RealVect xiRotated = rotateFwd(a_Xi);

  RealVect x = xiRotated;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (m_blockIndex[idir] == 1)
        {
          x[idir] -= 1.;
        }
    }
  return x;
}

RealVect
CartesianBlockRotateCS::mappedCoord(const RealVect& a_x) const
{
  RealVect xi = a_x;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      switch (m_blockIndex[idir])
        {
        case 0:
          if (a_x[idir] > 1.5)
            xi[idir] -= 2.;
          break;
        case 1:
          if (a_x[idir] < 0.5)
            xi[idir] += 3.;
          else
            xi[idir] += 1.;
          break;
        }
    }

  // Shift so that the block is centered at the origin,
  // then do the reverse rotation,
  // then shift back.
  RealVect xiCentered = xi - m_blockCenter;
  RealVect xiRotated;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      xiRotated[m_perm[idir]] = m_sign[idir] * xi[idir];
    }
  xiRotated += m_blockCenter;

  return xiRotated;
}

Real
CartesianBlockRotateCS::dXdXi(const RealVect& a_Xi,
                              int a_dirX,
                              int a_dirXi) const
{
  Real deriv = (a_dirX == m_perm[a_dirXi]) ? m_sign[a_dirXi] * 1. : 0.;
  return deriv;
}

void
CartesianBlockRotateCS::volFlux(FluxBox& a_volFlux,
                                const FluxBox& a_Nt,
                                const Box& a_box) const
{
  // This is adapted from CubedSphere2DPanelCS::volFlux().

  Box bx1 = grow(a_box, 1);
  CH_assert(a_Nt.box().contains(bx1));
  FluxBox NtX(bx1, SpaceDim);

  // When m_blockIndex[idir] == 0, signs[idir] == 1;
  // when m_blockIndex[idir] == 1, signs[idir] == -1.
  IntVect signs = IntVect::Unit - 2 * m_blockIndex;

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
          RealVect mappedLoc = rotateFwd(m_dx*iv + offset);
          for (int comp=m_volInterval.begin(); comp<=m_volInterval.end(); comp++)
            {
              // only use components in volInterval
              thisNtXdir(iv, comp) = signs[comp] *
                (mappedLoc[comp] - m_blockCenter[comp]) / Real(m_volInterval.size());
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

RealVect
CartesianBlockRotateCS::rotateFwd(const RealVect& a_xi) const
{
  // Shift so that the block is centered at the origin,
  // then do the rotation,
  // then shift back.
  RealVect xiCentered = a_xi - m_blockCenter;
  RealVect xiRotated;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      xiRotated[idir] = m_sign[idir] * xiCentered[m_perm[idir]];
    }
  xiRotated += m_blockCenter;
  return xiRotated;
}

RealVect
CartesianBlockRotateCS::rotateBack(const RealVect& a_xi) const
{
  // Shift so that the block is centered at the origin,
  // then do the reverse rotation,
  // then shift back.
  RealVect xiCentered = a_xi - m_blockCenter;
  RealVect xiRotated;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      xiRotated[m_perm[idir]] = m_sign[idir] * xiCentered[idir];
    }
  xiRotated += m_blockCenter;
  return xiRotated;
}

#include "NamespaceFooter.H"
