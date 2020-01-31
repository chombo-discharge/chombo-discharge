#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CubedSphere2DPanelCS.H"
#include "BoxIterator.H"
#include <cmath>

#include "CubedSphere2DF_F.H"
#include "FourthOrderUtil.H"
#include "FourthOrderUtilF_F.H"
#include "FCDivergenceF_F.H"

#include "NamespaceHeader.H"
// #include "UsingNamespace.H"

CubedSphere2DPanelCS::CubedSphere2DPanelCS(
    int nPanel,
    RealVect& dx,
    IntVect& ix
)
{
  CH_assert((nPanel >= 0) && (nPanel < 6));
  m_nPanel = nPanel;
  m_dx = dx;
  m_ix = ix;
  m_realDim = 3;
  m_flatMap = false;
}

CubedSphere2DPanelCS::~CubedSphere2DPanelCS()
{
}

void
CubedSphere2DPanelCS::setFlatMap(bool a_flatMap)
{
  m_flatMap = a_flatMap;
  if (m_flatMap)
    if (SpaceDim == 2)
      m_realDim = 2;
}

RealVect
CubedSphere2DPanelCS::realCoord(const RealVect& a_Xi) const
{
#if CH_SPACEDIM == 2
  // Real coordinate cannot be stored in a RealVect
  CH_assert(false);
#elif CH_SPACEDIM == 3
  Real vecXYZ[3];
  pointTransformEquiangularToCartesian(a_Xi, vecXYZ);
  RealVect XYZ;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      XYZ[idir] = vecXYZ[idir];
    }
  return XYZ;
#endif
/*
  // Take advantage of periodicity of tan function

  // Gnomonic coordinate equivalents
  Real dGX = tan(a_Xi[0]);
  Real dGY = tan(a_Xi[1]);

  // Delta parameter
  Real dInvDelta = 1.0 / sqrt(1.0 + dGX * dGX + dGY * dGY);

  RealVect xyzLoc;
  switch (m_nPanel)
  {
    case 0:
      xyzLoc[0] = dInvDelta;
      xyzLoc[1] = dGX * dInvDelta;
      xyzLoc[2] = dGY * dInvDelta;
      break;

    case 1:
      xyzLoc[0] = - dGX * dInvDelta;
      xyzLoc[1] = dInvDelta;
      xyzLoc[2] = dGY * dInvDelta;
      break;

    case 2:
      xyzLoc[0] = - dInvDelta;
      xyzLoc[1] = - dGX * dInvDelta;
      xyzLoc[2] = dGY * dInvDelta;
      break;

    case 3:
      xyzLoc[0] = dGX * dInvDelta;
      xyzLoc[1] = - dInvDelta;
      xyzLoc[2] = dGY * dInvDelta;
      break;

    case 4:
      xyzLoc[0] = - dGY * dInvDelta;
      xyzLoc[1] = dGX * dInvDelta;
      xyzLoc[2] = dInvDelta;
      break;

    case 5:
      xyzLoc[0] = dGY * dInvDelta;
      xyzLoc[1] = dGX * dInvDelta;
      xyzLoc[2] = - dInvDelta;
      break;

  }

  return xyzLoc;
*/
  return RealVect::Zero;
}

// given coordinates in mapped space, return locations in real space
void
CubedSphere2DPanelCS::realCoord(FArrayBox& a_x,
                                const FArrayBox& a_Xi,
                                const Box& a_box) const
{
  if (m_flatMap)
    {
      // flat view

      // I hate this kludge, but it's the only way I know how to get
      // the correct visualization.
      const IntVect& ivLo = a_box.smallEnd();
      const IntVect& ivHi = a_box.bigEnd();
      IntVect ivMid = (ivLo + ivHi)/2;
      Real midLon;
      {
        IntVect iv = ivMid;

        RealVect vecXi;
        vecXi[0] = a_Xi(iv, 0);
        vecXi[1] = a_Xi(iv, 1);
        RealVect rllXi;
        pointTransformEquiangularToLonLat(vecXi, rllXi);
        Real lon = rllXi[0];
        // Real lat = rllXi[1];
        if (lon < 0.) lon += 2 * M_PI;

        midLon = lon;
      }

      Real small = 0.001 * m_dx[0];
      // default implementation is probably inefficient, but it should work
      BoxIterator bit(a_box);
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();

          RealVect vecXi;
          vecXi[0] = a_Xi(iv, 0);
          vecXi[1] = a_Xi(iv, 1);
          RealVect rllXi;
          pointTransformEquiangularToLonLat(vecXi, rllXi);

          // Kludge for north pole.
          if (m_nPanel == 4)
            if (fabs(rllXi[1] - M_PI / 2.0) < small)
              {
                IntVect ivOff = iv - ivMid;
                // Set vecXi to be pulled back a small amount from the corner,
                // towards the middle of the box.
                for (int idir = 0; idir < SpaceDim; idir++)
                  {
                    vecXi[idir] -= (ivOff[idir] > 0) ? small : -small;
                  }
                pointTransformEquiangularToLonLat(vecXi, rllXi);
              }

          // Kludge for south pole.
          if (m_nPanel == 5)
            if (fabs(rllXi[1] + M_PI / 2.0) < small)
              {
                IntVect ivOff = iv - ivMid;
                // Set vecXi to be pulled back a small amount from the corner,
                // towards the middle of the box.
                for (int idir = 0; idir < SpaceDim; idir++)
                  {
                    vecXi[idir] -= (ivOff[idir] > 0) ? small : -small;
                  }
                pointTransformEquiangularToLonLat(vecXi, rllXi);
              }

          Real lon = rllXi[0];
          Real lat = rllXi[1];
          if (lon < 0.) lon += 2 * M_PI;

          // Kludge to keep longitude continuous within the box.
          if (lon < midLon - M_PI)
            lon += 2*M_PI;
          else if (lon > midLon + M_PI)
            lon -= 2*M_PI;

          a_x(iv, 0) = lon;
          a_x(iv, 1) = lat;
#if CH_SPACEDIM == 3
          Real rreal;
          Real rMap = a_Xi(iv, 2);
          FORT_SUBRMAPPEDTOREAL(CHF_REAL(rreal),
                                CHF_CONST_REAL(rMap));
          a_x(iv, 2) = rreal;
#endif
        }
    }
  else // !m_flatMap
    {
      FORT_CUBEDSPHERE2DMAPPEDTOREAL(CHF_FRA(a_x),
                                     CHF_CONST_FRA(a_Xi),
                                     CHF_CONST_INT(m_nPanel),
                                     CHF_BOX(a_box),
                                     CHF_CONST_REALVECT(m_dx),
                                     CHF_CONST_INTVECT(m_ix));
    }
}


// given coordinate in real space, return its location in the mapped space
RealVect
CubedSphere2DPanelCS::mappedCoord(const RealVect& a_x) const
{
  // Real coordinate cannot be stored in a RealVect
  CH_assert(false);
/*
  RealVect mappedXi;

  switch (m_nPanel)
  {
    case 0:
      mappedXi[0] = atan(a_x[1] / a_x[0]);
      mappedXi[1] = atan(a_x[2] / a_x[0]);
      break;

    case 1:
      mappedXi[0] = - atan(a_x[0] / a_x[1]) + M_PI;
      mappedXi[1] = atan(a_x[2] / a_x[1]);
      break;

    case 2:
      mappedXi[0] = atan(a_x[1] / a_x[0]) + 2.0 * M_PI;
      mappedXi[1] = - atan(a_x[2] / a_x[0]);
      break;

    case 3:
      mappedXi[0] = - atan(a_x[0] / a_x[1]) + 3.0 * M_PI;
      mappedXi[1] = - atan(a_x[2] / a_x[1]);
      break;

    case 4:
      mappedXi[0] = atan(a_x[1] / a_x[2]) + 4.0 * M_PI;
      mappedXi[1] = - atan(a_x[0] / a_x[2]);
      break;

    case 5:
      mappedXi[0] = - atan(a_x[1] / a_x[2]) + 5.0 * M_PI;
      mappedXi[1] = - atan(a_x[0] / a_x[2]);
      break;
  }

  return mappedXi;
*/
  return RealVect::Zero;
}

Real
CubedSphere2DPanelCS::dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const
{
  Real value = 0.0;

  // Gnomonic coordinates
  //RealVect a_Xi = mappedCoord(a_X);
  Real dGX = tan(a_Xi[0] - 0.25 * M_PI);
  Real dGY = tan(a_Xi[1] - 0.25 * M_PI);

  Real dInvDelta = 1.0 / sqrt(1.0 + dGX * dGX + dGY * dGY);
  Real dInvDelta3 = dInvDelta * dInvDelta * dInvDelta;

  if (a_dirX == 0)
    {
      // dx/dalpha
      if (a_dirXi == 0)
        {
          switch (m_nPanel)
            {
            case 0:
              value = (1.0 + dGX * dGX) * (- dGX * dInvDelta3);
              break;
            case 1:
              value = (1.0 + dGX * dGX) * (- (1.0 + dGY * dGY) * dInvDelta3);
              break;
            case 2:
              value = (1.0 + dGX * dGX) * (dGX * dInvDelta3);
              break;
            case 3:
              value = (1.0 + dGX * dGX) * ((1.0 + dGY * dGY) * dInvDelta3);
              break;
            case 4:
              value = (1.0 + dGX * dGX) * (dGX * dGY * dInvDelta3);
              break;
            case 5:
              value = (1.0 + dGX * dGX) * (- dGX * dGY * dInvDelta3);
              break;
            }
        }
      // dx/dbeta
      else if (a_dirXi == 1)
        {
          switch (m_nPanel)
            {
            case 0:
              value = (1.0 + dGY * dGY) * (- dGY * dInvDelta3);
              break;
            case 1:
              value = (1.0 + dGY * dGY) * (dGX * dGY * dInvDelta3);
              break;
            case 2:
              value = (1.0 + dGY * dGY) * (dGY * dInvDelta3);
              break;
            case 3:
              value = (1.0 + dGY * dGY) * (- dGX * dGY * dInvDelta3);
              break;
            case 4:
              value = (1.0 + dGY * dGY) * (- (1.0 + dGX * dGX) * dInvDelta3);
              break;
            case 5:
              value = (1.0 + dGY * dGY) * ((1.0 + dGX * dGX) * dInvDelta3);
              break;
            }
        }
    }
  else if (a_dirX == 1)
    {
      // dy/dalpha
      if (a_dirXi == 0)
        {
          switch (m_nPanel)
            {
            case 0:
              value = (1.0 + dGX * dGX) * ((1.0 + dGY * dGY) * dInvDelta3);
              break;
            case 1:
              value = (1.0 + dGX * dGX) * (- dGX * dInvDelta3);
              break;
            case 2:
              value = (1.0 + dGX * dGX) * (- (1.0 + dGY * dGY) * dInvDelta3);
              break;
            case 3:
              value = (1.0 + dGX * dGX) * (dGX * dInvDelta3);
              break;
            case 4:
              value = (1.0 + dGX * dGX) * ((1.0 + dGY * dGY) * dInvDelta3);
              break;
            case 5:
              value = (1.0 + dGX * dGX) * ((1.0 + dGY * dGY) * dInvDelta3);
              break;
            }
        }
      // dy/dbeta
      else if (a_dirXi == 1)
        {
          switch (m_nPanel)
            {
            case 0:
              value = (1.0 + dGY * dGY) * (- dGX * dGY * dInvDelta3);
              break;
            case 1:
              value = (1.0 + dGY * dGY) * (- dGY * dInvDelta3);
              break;
            case 2:
              value = (1.0 + dGY * dGY) * (dGX * dGY * dInvDelta3);
              break;
            case 3:
              value = (1.0 + dGY * dGY) * (dGY * dInvDelta3);
              break;
            case 4:
              value = (1.0 + dGY * dGY) * (- dGX * dGY * dInvDelta3);
              break;
            case 5:
              value = (1.0 + dGY * dGY) * (- dGX * dGY * dInvDelta3);
              break;
            }
        }
        

    }
  else if (a_dirX == 2)
    {
      // dz/dalpha
      if (a_dirXi == 0)
        {
          switch (m_nPanel)
            {
            case 0:
            case 1:
            case 2:
            case 3:
              value = (1.0 + dGX * dGX) * (- dGX * dGY * dInvDelta3);
              break;
            case 4:
              value = (1.0 + dGX * dGX) * (- dGX * dInvDelta3);
              break;
            case 5:
              value = (1.0 + dGX * dGX) * (dGX * dInvDelta3);
              break;
            }
        }
      // dz/dbeta
      else if (a_dirXi == 1)
        {
          switch (m_nPanel)
            {
            case 0:
            case 1:
            case 2:
            case 3:
              value = (1.0 + dGY * dGY) * ((1.0 + dGX * dGX) * dInvDelta3);
              break;
            case 4:
              value = (1.0 + dGY * dGY) * (- dGY * dInvDelta3);
              break;
            case 5:
              value = (1.0 + dGY * dGY) * (dGY * dInvDelta3);
              break;
            }
        }
    }
#if CH_SPACEDIM == 3
  Real rReal;
  Real rMap = a_Xi[2];
  FORT_SUBRMAPPEDTOREAL(CHF_REAL(rReal),
                        CHF_CONST_REAL(rMap));

  if (a_dirXi == 2)
    {
      RealVect xyz = realCoord(a_Xi);
      Real derivr;
      FORT_SUBDERIVRMAPPEDTOREAL(CHF_REAL(derivr),
                                 CHF_CONST_REAL(rMap));

      value = (xyz[a_dirX] / rReal) * derivr;
    }
  else
    { // Already done above, but assumed radius of 1.
      value *= rReal;
    }
#endif
  return value;
}

void
CubedSphere2DPanelCS::getNodeRealCoordinates(
                                             FArrayBox& a_nodeCoords,
                                             const Box& a_box
                                             ) const
{
  FORT_CUBEDSPHERE2DNODEREAL(CHF_FRA(a_nodeCoords),
                             CHF_CONST_INT(m_nPanel),
                             CHF_BOX(a_box),
                             CHF_CONST_REALVECT(m_dx),
                             CHF_CONST_INTVECT(m_ix));
}


void 
CubedSphere2DPanelCS::getFaceMappedCoordinates(
    FArrayBox& a_faceCoords,
    const int a_dir,
    const Box& a_box) const
{
  FORT_CUBEDSPHERE2DFACEMAPPED(CHF_FRA(a_faceCoords),
                               CHF_CONST_INT(m_nPanel),
                               CHF_BOX(a_box),
                               CHF_CONST_REALVECT(m_dx),
                               CHF_CONST_INT(a_dir),
                               CHF_CONST_INTVECT(m_ix));
}


void CubedSphere2DPanelCS::getCellMappedCoordinates(
                                                    FArrayBox& a_cellCoords,
                                                    const Box& a_box
                                                    ) const
{
  FORT_CUBEDSPHERE2DCELLMAPPED(CHF_FRA(a_cellCoords),
                               CHF_CONST_INT(m_nPanel),
                               CHF_BOX(a_box),
                               CHF_CONST_REALVECT(m_dx),
                               CHF_CONST_INTVECT(m_ix));
}

void CubedSphere2DPanelCS::getNodeMappedCoordinates(
                                                    FArrayBox& a_nodeCoords,
                                                    const Box& a_box
                                                    ) const
{
  FORT_CUBEDSPHERE2DNODEMAPPED(CHF_FRA(a_nodeCoords),
                               CHF_CONST_INT(m_nPanel),
                               CHF_BOX(a_box),
                               CHF_CONST_REALVECT(m_dx),
                               CHF_CONST_INTVECT(m_ix));
}


//-----------------------------------------------------------------------
void
CubedSphere2DPanelCS::
volFlux(FluxBox& a_volFlux,
        const FluxBox& a_Nt,
        const Box& a_box) const
{
  // Note that X needs to have one more ghost cell
  Box bx1 = grow(a_box, 1);
  CH_assert(a_Nt.box().contains(bx1));
  FluxBox NtX(bx1, SpaceDim);

  // Are these ever anything other than 0 and 1?
  // int compVol0 = m_volInterval.begin();
  // int compVol1 = m_volInterval.end();

  int signX = 1;
  switch (panel())
    {
    case 1:
    case 3:
    case 5:
      signX = -1;
    }

  int signY = 1;
  if (panel() > 1)
    {
      signY = -1;
    }

  // Compute the point values of NtX in the (alpha, beta) basis.
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
        mappedLoc += offset; // mapped coordinates of face center

        // int signX = (sin(mappedLoc[0] + 0.25 * M_PI) < 0.) ? -1 : 1;
        Real dGX = signX*tan(mappedLoc[0] - 0.25 * M_PI);
        // int signY = ((mappedLoc[0] + 0.25 * M_PI) > (2 * M_PI)) ? -1 : 1;
        Real dGY = signY*tan(mappedLoc[1] - 0.25 * M_PI);
        Real dDelta = sqrt(1.0 + dGX*dGX + dGY*dGY);
        thisNtXdir(iv, 0) = 0.5 * dGX / dDelta;
        thisNtXdir(iv, 1) = 0.5 * dGY / dDelta;
#if CH_SPACEDIM == 3
        Real rMap = mappedLoc[2];
        Real rReal;
        FORT_SUBRMAPPEDTOREAL(CHF_REAL(rReal),
                              CHF_CONST_REAL(rMap));
        thisNtXdir(iv, 0) *= rReal;
        thisNtXdir(iv, 1) *= rReal;
        thisNtXdir(iv, 2) = rReal * 0.5 / dDelta;
#endif
      }
  } // end loop over directions

  // Convert point values to 4th-order face averages
  fourthOrderAverageFace(NtX);
  // added by petermc, 21 Nov 2011:  a_volFlux[dir] = <X[dir]>
  for (int dir=0; dir<SpaceDim; dir++)
  {
    // a_volFlux[dir][0] := NtX[dir][dir]
    a_volFlux[dir].copy(NtX[dir], dir, 0);
  }

  // removed by petermc, 21 Nov 2011
  // Compute the volume flux on the faces.
  //  computeMetricTermProductAverage(a_volFlux, X, a_Nt, SpaceDim, X, a_box);
}
//-----------------------------------------------------------------------

void
CubedSphere2DPanelCS::cellVol(FArrayBox& a_vol,
                              const FluxBox& a_N,
                              const Box& a_box) const
{
  FORT_CUBEDSPHERE2DCELLVOL(CHF_FRA1(a_vol, 0),
                            CHF_CONST_INT(m_nPanel),
                            CHF_BOX(a_box),
                            CHF_CONST_REALVECT(m_dx),
                            CHF_CONST_INTVECT(m_ix));
}

void CubedSphere2DPanelCS::getN(FluxBox& a_N, const Box& a_box) const
{
  a_N.setVal(0.0);
  // a_N[0] and a_N[1] both have 2*2 = 4 components
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      FArrayBox& Nfab = a_N[idir];
      // added by petermc, 12 Apr 2011
      Box faceBox(a_box);
      faceBox.surroundingNodes(idir);
      FORT_CUBEDSPHERE2DNORMAL(CHF_FRA(Nfab),
                               CHF_CONST_INT(m_nPanel),
                               CHF_BOX(faceBox),
                               CHF_CONST_INT(idir),
                               CHF_CONST_REALVECT(m_dx),
                               CHF_CONST_INTVECT(m_ix));
    }
}

/// computes cell-averaged J
void CubedSphere2DPanelCS::getAvgJ(FArrayBox& a_avgJ,
                                   const Box& a_box) const
{
  // we don't need N the way we're computing this, but it's required by
  // the API
  FluxBox Nbogus;

  cellVol(a_avgJ, Nbogus, a_box);

  Real mappedVol = m_dx.product();
  a_avgJ /= mappedVol;
}

/// computes cell-averaged J
void CubedSphere2DPanelCS::getAvgJ(FArrayBox& a_avgJ,
                                   const FluxBox& a_volFlux,
                                   const Box& a_box) const
{
  RealVect panelDxSign = m_dx;
  // Removed by petermc, 21 Nov 2011
  //  switch (panel())
  //    {
  //    case 1:
  //    case 3:
  //    case 5:
  //      panelDxSign[0] *= -1;
  //    }
  //  if (panel() > 1)
  //    {
  //      panelDxSign[1] *= -1;
  //    }

  // Set a_avgJ = div(a_volFlux).
  a_avgJ.setVal(0.0, a_box, 0, a_avgJ.nComp());
  for (int idir = 0; idir < SpaceDim; idir++)
    { // Computes a_avgJ += diff(a_volFlux[idir]) / panelDxSign[idir].
      FORT_FCDIVERGENCE(CHF_CONST_FRA(a_volFlux[idir]),
                        CHF_FRA(a_avgJ),
                        CHF_BOX(a_box),
                        CHF_CONST_REAL(panelDxSign[idir]),
                        CHF_INT(idir));
    }
}

/// computes cell-averaged 1/J
void CubedSphere2DPanelCS::getAvgJinverse(
                                          FluxBox& a_avgJinverse,
                                          const Box& a_box
                                          ) const
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      FArrayBox& avgJinvFab = a_avgJinverse[idir];
      FORT_CUBEDSPHERE2DAVGJINV(CHF_FRA1(avgJinvFab, 0),
                                CHF_CONST_INT(m_nPanel),
                                CHF_BOX(a_box),
                                CHF_CONST_INT(idir),
                                CHF_CONST_REALVECT(m_dx),
                                CHF_CONST_INTVECT(m_ix));
    }
}

/// Jacobian evaluated at location X in mapped space
Real CubedSphere2DPanelCS::pointwiseJ(const RealVect& a_Xi) const
{
  Real dGX = tan(a_Xi[0] - 0.25 * M_PI);
  Real dGY = tan(a_Xi[1] - 0.25 * M_PI);

  Real dDelta = sqrt(1.0 + dGX * dGX + dGY * dGY);

  Real J = (1.0 + dGX * dGX) * (1.0 + dGY * dGY) / (dDelta * dDelta * dDelta);
#if CH_SPACEDIM == 3
  Real rMap = a_Xi[2];
  Real derivr;
  FORT_SUBDERIVRMAPPEDTOREAL(CHF_REAL(derivr),
                             CHF_CONST_REAL(rMap));
  J *= derivr;
#endif  
  return J;
}

/// Jacobian evaluated at locations Xi in mapped space
void CubedSphere2DPanelCS::pointwiseJ(
                                      FArrayBox& a_J,
                                      const FArrayBox& a_Xi,
                                      const Box& a_box
                                      ) const
{
  FORT_CUBEDSPHERE2DPOINTJ(CHF_FRA1(a_J, 0),
                           CHF_CONST_FRA(a_Xi),
                           CHF_CONST_INT(m_nPanel),
                           CHF_BOX(a_box),
                           CHF_CONST_REALVECT(m_dx),
                           CHF_CONST_INTVECT(m_ix));
}

/// returns integral of divergence over mapped-grid cells
void CubedSphere2DPanelCS::computeDivergence(
                                             FArrayBox& a_divF,
                                             const FluxBox& a_F,
                                             const FluxBox& a_N,
                                             const Box& a_box,
                                             Interval& divInterval
) const
{
}

/// transform a point from mapped-coordinate basis to latitude-longitude
/// coordinate basis
void CubedSphere2DPanelCS::pointTransformEquiangularToLonLat(
  const RealVect& a_xi,
  RealVect& a_rllXi
) const
{
  Real eps = DBL_EPSILON;
  // In 3D, component 2 is just copied from a_xi to a_rllXi.
  FORT_CUBEDSPHERE2DEQUIANGULARTOLONLAT(CHF_REALVECT(a_rllXi),
                                        CHF_CONST_REALVECT(a_xi),
                                        CHF_CONST_INT(m_nPanel),
                                        CHF_CONST_REAL(eps));
}


/// transform a point from mapped-coordinate basis to Cartesian
/// coordinate basis
void
CubedSphere2DPanelCS::pointTransformEquiangularToCartesian(
                                                           const RealVect& a_xi,
                                                           Real * a_xyz)
  const
{
  FORT_CUBEDSPHERE2DEQUIANGULARTOCARTESIAN(CHF_R1D(a_xyz, 3),
                                           CHF_CONST_REALVECT(a_xi),
                                           CHF_CONST_INT(m_nPanel));
}


/// transform a FAB of points from mapped-coordinate basis to latitude-longitude
/// coordinate basis
void CubedSphere2DPanelCS::fabTransformEquiangularToLonLat(
                                                           const FArrayBox& a_xiFab,
                                                           FArrayBox& a_rllXiFab
                                                           ) const
{
  CH_assert(a_xiFab.box() == a_rllXiFab.box());
  Real eps = DBL_EPSILON;
  // In 3D, component 2 is just copied from a_xiFab to a_rllXiFab.
  FORT_CUBEDSPHERE2DFABEQUIANGULARTOLONLAT(CHF_FRA(a_rllXiFab),
                                           CHF_CONST_FRA(a_xiFab),
                                           CHF_CONST_INT(m_nPanel),
                                           CHF_CONST_REAL(eps));
}


/// transform a FAB of points from mapped-coordinate basis to Cartesian
/// coordinate basis
void
CubedSphere2DPanelCS::fabTransformEquiangularToCartesian(
                                                         const FArrayBox& a_xiFab,
                                                         FArrayBox& a_xyzFab) const
{
  CH_assert(a_xiFab.box() == a_xyzFab.box());
  FORT_CUBEDSPHERE2DFABEQUIANGULARTOCARTESIAN(CHF_FRA(a_xyzFab),
                                              CHF_CONST_FRA(a_xiFab),
                                              CHF_CONST_INT(m_nPanel));
}


/// transform a vector from mapped-coordinate basis to real-coordinate basis
void CubedSphere2DPanelCS::vectorTransformEquiangularToCartesian(
  const RealVect& a_xi,
  const Real * a_vecCS,
  Real * a_vecXYZ
) const
{
  CH_assert((m_nPanel >= 0) && (m_nPanel <= 5));

  // Gnomonic coordinate equivalents
  Real dGX = tan(a_xi[0] - 0.25 * M_PI);
  Real dGY = tan(a_xi[1] - 0.25 * M_PI);

  // Delta parameter
  Real dInvDelta  = 1.0 / sqrt(1.0 + dGX * dGX + dGY * dGY);
  Real dInvDelta3 = dInvDelta * dInvDelta * dInvDelta;

  switch (m_nPanel)
  {
    case 0:
      a_vecXYZ[0] =
        - dGX * (1.0 + dGX * dGX) * dInvDelta3 * a_vecCS[0]
        - dGY * (1.0 + dGY * dGY) * dInvDelta3 * a_vecCS[1];

      a_vecXYZ[1] =
          (1.0 + dGX * dGX) * (1.0 + dGY * dGY) * dInvDelta3 * a_vecCS[0]
        - (1.0 + dGY * dGY) * dGX * dGY * dInvDelta3 * a_vecCS[1];

      a_vecXYZ[2] =
        - (1.0 + dGX * dGX) * dGX * dGY * dInvDelta3 * a_vecCS[0]
        + (1.0 + dGX * dGX) * (1.0 + dGY * dGY) * dInvDelta3 * a_vecCS[1];
      break;

    case 1:
      a_vecXYZ[0] =
        - (1.0 + dGX * dGX) * (1.0 + dGY * dGY) * dInvDelta3 * a_vecCS[0]
        + (1.0 + dGY * dGY) * dGX * dGY * dInvDelta3 * a_vecCS[1];

      a_vecXYZ[1] =
        - dGX * (1.0 + dGX * dGX) * dInvDelta3 * a_vecCS[0]
        - dGY * (1.0 + dGY * dGY) * dInvDelta3 * a_vecCS[1];

      a_vecXYZ[2] =
        - (1.0 + dGX * dGX) * dGX * dGY * dInvDelta3 * a_vecCS[0]
        + (1.0 + dGX * dGX) * (1.0 + dGY * dGY) * dInvDelta3 * a_vecCS[1];
      break;

    case 2:
      a_vecXYZ[0] =
          dGX * (1.0 + dGX * dGX) * dInvDelta3 * a_vecCS[0]
        + dGY * (1.0 + dGY * dGY) * dInvDelta3 * a_vecCS[1];

      a_vecXYZ[1] =
        - (1.0 + dGX * dGX) * (1.0 + dGY * dGY) * dInvDelta3 * a_vecCS[0]
        + (1.0 + dGY * dGY) * dGX * dGY * dInvDelta3 * a_vecCS[1];

      a_vecXYZ[2] =
        - (1.0 + dGX * dGX) * dGX * dGY * dInvDelta3 * a_vecCS[0]
        + (1.0 + dGX * dGX) * (1.0 + dGY * dGY) * dInvDelta3 * a_vecCS[1];
      break;

    case 3:
      a_vecXYZ[0] =
          (1.0 + dGX * dGX) * (1.0 + dGY * dGY) * dInvDelta3 * a_vecCS[0]
        - (1.0 + dGY * dGY) * dGX * dGY * dInvDelta3 * a_vecCS[1];

      a_vecXYZ[1] =
          dGX * (1.0 + dGX * dGX) * dInvDelta3 * a_vecCS[0]
        + dGY * (1.0 + dGY * dGY) * dInvDelta3 * a_vecCS[1];

      a_vecXYZ[2] =
        - dGX * dGY * (1.0 + dGX * dGX) * dInvDelta3 * a_vecCS[0]
        + (1.0 + dGX * dGX) * (1.0 + dGY * dGY) * dInvDelta3 * a_vecCS[1];
      break;

    case 4:
      a_vecXYZ[0] =
          (1.0 + dGX * dGX) * dGX * dGY * dInvDelta3 * a_vecCS[0]
        - (1.0 + dGX * dGX) * (1.0 + dGY * dGY) * dInvDelta3 * a_vecCS[1];

      a_vecXYZ[1] =
          (1.0 + dGX * dGX) * (1.0 + dGY * dGY) * dInvDelta3 * a_vecCS[0]
        - (1.0 + dGY * dGY) * dGX * dGY * dInvDelta3 * a_vecCS[1];

      a_vecXYZ[2] =
        - dGX * (1.0 + dGX * dGX) * dInvDelta3 * a_vecCS[0]
        - dGY * (1.0 + dGY * dGY) * dInvDelta3 * a_vecCS[1];
      break;

    case 5:
      a_vecXYZ[0] =
        - dGX * dGY * (1.0 + dGX * dGX) * dInvDelta3 * a_vecCS[0]
        + (1.0 + dGX * dGX) * (1.0 + dGY * dGY) * dInvDelta3 * a_vecCS[1];

      a_vecXYZ[1] =
          (1.0 + dGX * dGX) * (1.0 + dGY * dGY) * dInvDelta3 * a_vecCS[0]
        - (1.0 + dGY * dGY) * dGX * dGY * dInvDelta3 * a_vecCS[1];

      a_vecXYZ[2] =
          dGX * (1.0 + dGX * dGX) * dInvDelta3 * a_vecCS[0]
        + dGY * (1.0 + dGY * dGY) * dInvDelta3 * a_vecCS[1];
      break;
  }
#if CH_SPACEDIM == 3
  // Up to now, we have ignored the component a_vecCS[2] * rhat.
  // Note rhat = xu*ihat + yu*jhat + zu*khat,
  // where xu = x/r, yu = y/r, zu = z/r.
  Real rMap = a_xi[2];
  Real rReal;
  FORT_SUBRMAPPEDTOREAL(CHF_REAL(rReal),
                        CHF_CONST_REAL(rMap));
  RealVect xyz = realCoord(a_xi);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_vecXYZ[idir] += (xyz[idir] / rReal) * a_vecCS[2];
    }
#endif
}

/// transform a vector from real-coordinate basis to mapped-coordinate basis
void CubedSphere2DPanelCS::vectorTransformCartesianToEquiangular(
  const RealVect& a_xi,
  const Real * a_vecXYZ,
  Real * a_vecCS
) const
{
  // I don't think this is called anymore.
  CH_assert(false);

  CH_assert((m_nPanel >= 0) && (m_nPanel <= 5));

  // Gnomonic coordinate equivalents
  Real dGX = tan(a_xi[0] - 0.25 * M_PI);
  Real dGY = tan(a_xi[1] - 0.25 * M_PI);

  // Delta parameter
  Real dDelta  = sqrt(1.0 + dGX * dGX + dGY * dGY);

  switch (m_nPanel)
  {
    case 0:
      a_vecCS[0] =
        - dGX / (1.0 + dGX * dGX) * dDelta * a_vecXYZ[0]
        + 1.0 / (1.0 + dGX * dGX) * dDelta * a_vecXYZ[1];

      a_vecCS[1] =
        - dGY / (1.0 + dGY * dGY) * dDelta * a_vecXYZ[0]
        + 1.0 / (1.0 + dGY * dGY) * dDelta * a_vecXYZ[2];
      break;

    case 1:
      a_vecCS[0] =
        - 1.0 / (1.0 + dGX * dGX) * dDelta * a_vecXYZ[0]
        - dGX / (1.0 + dGX * dGX) * dDelta * a_vecXYZ[1];

      a_vecCS[1] =
        - dGY / (1.0 + dGY * dGY) * dDelta * a_vecXYZ[1]
        + 1.0 / (1.0 + dGY * dGY) * dDelta * a_vecXYZ[2];
      break;

    case 2:
      a_vecCS[0] =
          dGX / (1.0 + dGX * dGX) * dDelta * a_vecXYZ[0]
        - 1.0 / (1.0 + dGX * dGX) * dDelta * a_vecXYZ[1];

      a_vecCS[1] =
          dGY / (1.0 + dGY * dGY) * dDelta * a_vecXYZ[0]
        + 1.0 / (1.0 + dGY * dGY) * dDelta * a_vecXYZ[2];
      break;

    case 3:
      a_vecCS[0] =
          1.0 / (1.0 + dGX * dGX) * dDelta * a_vecXYZ[0]
        + dGX / (1.0 + dGX * dGX) * dDelta * a_vecXYZ[1];

      a_vecCS[1] =
          dGY / (1.0 + dGY * dGY) * dDelta * a_vecXYZ[1]
        + 1.0 / (1.0 + dGY * dGY) * dDelta * a_vecXYZ[2];
      break;

    case 4:
      a_vecCS[0] =
          1.0 / (1.0 + dGX * dGX) * dDelta * a_vecXYZ[1]
        - dGX / (1.0 + dGX * dGX) * dDelta * a_vecXYZ[2];

      a_vecCS[1] =
        - 1.0 / (1.0 + dGY * dGY) * dDelta * a_vecXYZ[0]
        - dGY / (1.0 + dGY * dGY) * dDelta * a_vecXYZ[2];
      break;

    case 5:
      a_vecCS[0] =
          1.0 / (1.0 + dGX * dGX) * dDelta * a_vecXYZ[1]
        + dGX / (1.0 + dGX * dGX) * dDelta * a_vecXYZ[2];

      a_vecCS[1] =
          1.0 / (1.0 + dGY * dGY) * dDelta * a_vecXYZ[0]
        + dGY / (1.0 + dGY * dGY) * dDelta * a_vecXYZ[2];
      break;
  }
}

/// transform a vector from latitude-longitude coordinate basis to
/// mapped-coordinate basis
void CubedSphere2DPanelCS::vectorTransformLatLonToEquiangular(
  const RealVect& a_xi,
  const Real * a_vecRLL,
  Real * a_vecCS
) const
{
  CH_assert((m_nPanel >= 0) && (m_nPanel <= 5));

  // Gnomonic coordinate equivalents
  Real dGX = tan(a_xi[0] - 0.25 * M_PI);
  Real dGY = tan(a_xi[1] - 0.25 * M_PI);

  // Delta coordinate parameter
  Real dDelta2 = 1.0 + dGX * dGX + dGY * dGY;

  // Latitude of point
  Real lat;

  // Radius of point in gnomonic coordinates
  Real gnoRadius;

  // Latitudinal and longitudinal component of vector (scaled)
  Real vecLon;
  Real vecLat = a_vecRLL[1];

  switch(m_nPanel)
  {
    // Equatorial panels
    case 0:
    case 1:
    case 2:
    case 3:
      lat = atan(dGY / sqrt(1.0 + dGX * dGX));
      vecLon = a_vecRLL[0] / cos(lat);

      a_vecCS[0] = vecLon;
      a_vecCS[1] =
        dGX * dGY / (1.0 + dGY * dGY) * vecLon
        + dDelta2 / ((1.0 + dGY * dGY) * sqrt(1.0 + dGX * dGX)) * vecLat;
      break;

    // North polar panel
    case 4:
      lat = 0.5 * M_PI - atan(sqrt(dGX * dGX + dGY * dGY));
      vecLon = a_vecRLL[0] / cos(lat);

      gnoRadius = sqrt(dGX * dGX + dGY * dGY);

      a_vecCS[0] =
        - dGY / (1.0 + dGX * dGX) * vecLon
        - dDelta2 * dGX / ((1.0 + dGX * dGX) * gnoRadius) * vecLat;

      a_vecCS[1] =
        dGX / (1.0 + dGY * dGY) * vecLon
        - dDelta2 * dGY / ((1.0 + dGY * dGY) * gnoRadius) * vecLat;
      break;

     // South polar panel
     case 5:
       lat = -0.5 * M_PI + atan(sqrt(dGX * dGX + dGY * dGY));
       vecLon = a_vecRLL[0] / cos(lat);

       gnoRadius = sqrt(dGX * dGX + dGY * dGY);

      a_vecCS[0] =
        dGY / (1.0 + dGX * dGX) * vecLon
        + dDelta2 * dGX / ((1.0 + dGX * dGX) * gnoRadius) * vecLat;

      a_vecCS[1] =
        - dGX / (1.0 + dGY * dGY) * vecLon
        + dDelta2 * dGY / ((1.0 + dGY * dGY) * gnoRadius) * vecLat;
      break;
  }
#if CH_SPACEDIM == 3
  // Radial component does not change.
  a_vecCS[2] = a_vecRLL[2];
#endif
}

/// transform a FAB of vectors from latitude-longitude coordinate basis to
/// mapped-coordinate basis
void CubedSphere2DPanelCS::fabVectorTransformLatLonToEquiangular(
                                                              const FArrayBox& a_xiFab,
                                                              const FArrayBox& a_vecRLLFab,
                                                              FArrayBox& a_vecCSFab) const
{
  CH_assert((m_nPanel >= 0) && (m_nPanel <= 5));
  CH_assert(a_xiFab.box() == a_vecRLLFab.box());
  CH_assert(a_xiFab.box() == a_vecCSFab.box());
  // In 3D, component 2 is just copied from a_vecRLLFab to a_vecCSFab.
  FORT_CUBEDSPHERE2DFABVECTORLATLONTOEQUIANGULAR(CHF_CONST_FRA(a_xiFab),
                                                 CHF_CONST_FRA(a_vecRLLFab),
                                                 CHF_FRA(a_vecCSFab),
                                                 CHF_CONST_INT(m_nPanel));
}


/// transform a FAB of vectors from mapped-coordinate basis to
/// latitude-longitude coordinate basis
void CubedSphere2DPanelCS::fabVectorTransformEquiangularToLatLon(
                                                              const FArrayBox& a_xiFab,
                                                              const FArrayBox& a_vecCSFab,
                                                              FArrayBox& a_vecRLLFab) const
{
  CH_assert((m_nPanel >= 0) && (m_nPanel <= 5));
  CH_assert(a_xiFab.box() == a_vecRLLFab.box());
  CH_assert(a_xiFab.box() == a_vecCSFab.box());
  // In 3D, component 2 is just copied from a_vecCSFab to a_vecRLLFab.
  FORT_CUBEDSPHERE2DFABVECTOREQUIANGULARTOLATLON(CHF_CONST_FRA(a_xiFab),
                                                 CHF_CONST_FRA(a_vecCSFab),
                                                 CHF_FRA(a_vecRLLFab),
                                                 CHF_CONST_INT(m_nPanel));
}


/// transform a FAB of SpaceDim-vectors from mapped-coordinate basis to
/// real-coordinate basis at cell centers
void CubedSphere2DPanelCS::vectorTransformMappedToRealCenterFab(
                                                                FArrayBox& a_vectorFab
                                                                ) const
{
  Real vecCS[SpaceDim];
  Real vecXYZ[3];

  CH_assert(a_vectorFab.nComp() == 3); // use 2 comps on input, 3 on output
  BoxIterator bit(a_vectorFab.box());
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();

      //      RealVect ABcenter = m_dx * (iv - m_ix)
      //        + (dx() - 0.25 * M_PI) * RealVect::Unit;
      // changed by petermc, 12 Aug 2010
      // center of cell iv in mapped space
      RealVect ABcenter = iv * m_dx + 0.5 * m_dx;

      for (int idir = 0; idir < SpaceDim; idir++) // in 2D, no radial component
        {
          vecCS[idir] = a_vectorFab(iv, idir);
        }

      vectorTransformEquiangularToCartesian(ABcenter, vecCS, vecXYZ);

      for (int idir = 0; idir < 3; idir++) // always 3
        {
          a_vectorFab(iv, idir) = vecXYZ[idir];
        }
    }
}

/// transform a FAB of SpaceDim-vectors from real-coordinate basis to
/// mapped-coordinate basis at cell centers
void
CubedSphere2DPanelCS::vectorTransformRealToMappedCenterFab(
                                                           FArrayBox& a_vectorFab
) const
{
  CH_assert(a_vectorFab.nComp() == 3); // use 3 comps on input, 2 on output
  const Box& bx = a_vectorFab.box();
  FORT_CUBEDSPHERE2DVECTORREALTOMAPPEDCENTER(CHF_FRA(a_vectorFab),
                                             CHF_CONST_INT(m_nPanel),
                                             CHF_BOX(bx),
                                             CHF_CONST_REALVECT(m_dx),
                                             CHF_CONST_INTVECT(m_ix));
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
CubedSphere2DPanelCS::
getNMatrixEntry(const RealVect& a_Xi,
                int a_s, int a_d, int a_d1,
                int a_row, int a_column) const
{
  Real entry = 0.0;
  // simplest thing -- if row = s, then A_ij = delta_dj
  if (a_row == a_s)
    {
      if (a_column == a_d) entry = 1.0;
    }
  else if (a_column == a_d1)
    {
#if CH_SPACEDIM == 3
      // same as in NewFourthOrderCoordSys::getNMatrixEntry()
      entry = realCoord(a_Xi)[a_row];
#else
      // We can't use realCoord() here, since that method assumes that the
      // dimensions of the domain and the range spaces are the same.
      // Instead, we directly compute the (a_row)th component of X.
      RealVect ll;
      pointTransformEquiangularToLonLat(a_Xi, ll);
      entry = ll[a_row];
#endif
    }
  else
    {
      entry = dXdXi(a_Xi, a_row, a_column);
    }

  return entry;
}

//-----------------------------------------------------------------------
Real
CubedSphere2DPanelCS::getN(const RealVect& a_Xi,
                           int a_s, int a_d, int a_d1) const
{
#if CH_SPACEDIM == 3
  return NewFourthOrderCoordSys::getN(a_Xi, a_s, a_d, a_d1);
#else
  Real retval = 0.;
  if (a_s == a_d)
    {
      Real dGX = tan(a_Xi[0] - 0.25 * M_PI);
      Real dGY = tan(a_Xi[1] - 0.25 * M_PI);
      Real dDelta = sqrt(1.0 + dGX*dGX + dGY*dGY);
      if (a_d == 0)
        {
          retval = dGY / dDelta;
        }
      else if (a_d == 1)
        {
          retval = dGX / dDelta;
        }
      else
        {
          MayDay::Error("CubedSphere2DPanelCS::getN dimension must be 0 or 1");
        }
    }

  return retval;
#endif
}

//-----------------------------------------------------------------------

void
CubedSphere2DPanelCS::contravariantMetric(FArrayBox& a_metric,
                                          int a_dir) const
{
  CH_assert(a_metric.nComp() == SpaceDim);
  const Box& bx = a_metric.box();

  FArrayBox xiFab(bx, SpaceDim);
  getCenterMappedCoordinates(xiFab, bx); // according to centering of bx

  FORT_CUBEDSPHERE2DCONTRAVARIANTMETRIC(CHF_CONST_FRA(xiFab),
                                        CHF_FRA(a_metric),
                                        CHF_CONST_INT(a_dir));
}

//-----------------------------------------------------------------------

void
CubedSphere2DPanelCS::orthonormalize(
                                     const FArrayBox& a_csFab,
                                     FArrayBox& a_orthoFab,
                                     const Box& a_box,
                                     int a_idir,
                                     const IntVect& a_csComps,
                                     const IntVect& a_orthoComps) const
{
  FArrayBox xiFab(a_box, SpaceDim);
  getCenterMappedCoordinates(xiFab, a_box); // according to centering of a_box

  FORT_ORTHONORMALIZE(CHF_CONST_FRA(xiFab),
                      CHF_CONST_FRA(a_csFab),
                      CHF_FRA(a_orthoFab),
                      CHF_CONST_INT(a_idir),
                      CHF_CONST_INTVECT(a_csComps),
                      CHF_CONST_INTVECT(a_orthoComps));
}

//-----------------------------------------------------------------------

void
CubedSphere2DPanelCS::deorthonormalize(
                                       const FArrayBox& a_orthoFab,
                                       FArrayBox& a_csFab,
                                       const Box& a_box,
                                       int a_idir,
                                       const IntVect& a_orthoComps,
                                       const IntVect& a_csComps) const
{
  FArrayBox xiFab(a_box, SpaceDim);
  getCenterMappedCoordinates(xiFab, a_box); // according to centering of a_box

  FORT_DEORTHONORMALIZE(CHF_CONST_FRA(xiFab),
                        CHF_CONST_FRA(a_orthoFab),
                        CHF_FRA(a_csFab),
                        CHF_CONST_INT(a_idir),
                        CHF_CONST_INTVECT(a_orthoComps),
                        CHF_CONST_INTVECT(a_csComps));
}

//-----------------------------------------------------------------------

void
CubedSphere2DPanelCS::getOrthonormalizingMatrix(
                                                FArrayBox& a_matrixFab,
                                                const Box& a_box,
                                                int a_idir) const
{
  FArrayBox xiFab(a_box, SpaceDim);
  getCenterMappedCoordinates(xiFab, a_box); // according to centering of a_box

  FORT_ORTHONORMALIZINGMATRIX(CHF_CONST_FRA(xiFab),
                              CHF_FRA(a_matrixFab),
                              CHF_BOX(a_box),
                              CHF_CONST_INT(a_idir));
}

//-----------------------------------------------------------------------

void
CubedSphere2DPanelCS::getDeorthonormalizingMatrix(
                                                  FArrayBox& a_matrixFab,
                                                  const Box& a_box,
                                                  int a_idir) const
{
  FArrayBox xiFab(a_box, SpaceDim);
  getCenterMappedCoordinates(xiFab, a_box); // according to centering of a_box

  FORT_DEORTHONORMALIZINGMATRIX(CHF_CONST_FRA(xiFab),
                                CHF_FRA(a_matrixFab),
                                CHF_BOX(a_box),
                                CHF_CONST_INT(a_idir));
}

//-----------------------------------------------------------------------

void
CubedSphere2DPanelCS::curl(
                           const FArrayBox& a_vecFab,
                           FArrayBox& a_curlFab,
                           const Box& a_box) const
{
  Box bx1 = grow(a_box, 1);
  CH_assert(a_vecFab.box().contains(bx1));
  FArrayBox xiFab(bx1, SpaceDim);
  getCenterMappedCoordinates(xiFab, bx1); // according to centering of a_box

  // Get diffVecFab[idir] = d(a_vecFab)/dx[idir].
  // This will have SpaceDim components, because a_vecFab does.
  Tuple<FArrayBox, SpaceDim> diffVecFab;
  for (int idir=0; idir<SpaceDim; idir++)
    {
      diffVecFab[idir].define(a_box, SpaceDim);
      // FORT_UDIVCENTERGRAD returns undivided centered differences.
      FORT_UDIVCENTERGRAD(CHF_FRA(diffVecFab[idir]),
                          CHF_CONST_FRA(a_vecFab),
                          CHF_BOX(a_box),
                          CHF_CONST_INT(idir));
      diffVecFab[idir] /= m_dx[idir];
    }
  CH_assert(a_curlFab.box().contains(a_box));
  FORT_CUBEDSPHERE2DCURLR(CHF_CONST_FRA(xiFab),
                          CHF_CONST_FRA(diffVecFab[0]),
                          CHF_CONST_FRA(diffVecFab[1]),
                          CHF_FRA1(a_curlFab, 0),
                          CHF_BOX(a_box));
}

//-----------------------------------------------------------------------

void
CubedSphere2DPanelCS::curlSpherical(
                                    const FArrayBox& a_vecRLLFab,
                                    FArrayBox& a_curlFab,
                                    const Box& a_box) const
{
  Box bx1 = grow(a_box, 1);
  CH_assert(a_vecRLLFab.box().contains(bx1));
  FArrayBox xiFab(bx1, SpaceDim);
  getCenterMappedCoordinates(xiFab, bx1); // according to centering of a_box

  // int LON = 0;
  int LAT = 1;

  FArrayBox lonlatFab(bx1, SpaceDim);
  fabTransformEquiangularToLonLat(xiFab, lonlatFab);

  FArrayBox coslatFab(bx1, 1);
  BoxIterator bit1(bx1);
  for (bit1.begin(); bit1.ok(); ++bit1)
    {
      IntVect iv = bit1();
      coslatFab(iv, 0) = cos(lonlatFab(iv, LAT));
    }

  // Set vecModFab = (velZonal*cos(theta), velMeridional)
  // where theta is latitude = lonlatFab[1]
  int ZONALCOSLAT = 0;
  int MERIDIONAL = 1;
  FArrayBox vecModFab(bx1, SpaceDim);
  vecModFab.copy(a_vecRLLFab);
  vecModFab.mult(coslatFab, 0, ZONALCOSLAT);

  // Get diffVecModFab[idir] = d(vecModFab)/dx[idir].
  // This will have SpaceDim components, because vecModFab does.
  Tuple<FArrayBox, SpaceDim> diffVecModFab;
  for (int idir=0; idir<SpaceDim; idir++)
    {
      diffVecModFab[idir].define(a_box, SpaceDim);
      // FORT_UDIVCENTERGRAD returns undivided centered differences.
      FORT_UDIVCENTERGRAD(CHF_FRA(diffVecModFab[idir]),
                          CHF_CONST_FRA(vecModFab),
                          CHF_BOX(a_box),
                          CHF_CONST_INT(idir));
      diffVecModFab[idir] /= m_dx[idir];
    }
  // diffVecModFab[ALPHA] contains d(vecModFab)/dalpha
  // diffVecModFab[BETA] contains d(vecModFab)/dbeta
  // Each of these contains two components: ZONALCOSLAT and MERIDIONAL.

  // Now curl = ( d(vecModFab[1])/dlon - d(vecModFab[0])/dlat ) / cos(lat)

  // Chain rule:
  // dw/dlon = dw/dalpha * dalpha/dlon + dw/dbeta * dbeta/dlon
  // dw/dlat = dw/dalpha * dalpha/dlat + dw/dbeta * dbeta/dlat
  // where
  // diffVecModFab[0] contains dw/dalpha
  // diffVecModFab[1] contains dw/dbeta
  // so we need to find d{alpha, beta}/d{lon, lat}.
  // partialsFab contains:
  // dalpha/dlon, dbeta/dlon, dalpha/dlat, dbeta/dlat.
  FArrayBox partialsFab(a_box, SpaceDim*SpaceDim);
  FORT_CUBEDSPHERE2DPARTIALS(CHF_CONST_FRA(lonlatFab),
                             CHF_FRA(partialsFab),
                             CHF_CONST_INT(m_nPanel),
                             CHF_BOX(a_box));
  int ALPHALON = 0;
  int BETALON = 1;
  int ALPHALAT = 2;
  int BETALAT = 3;

  FArrayBox workFab(a_box, 1);

  int ALPHA = 0;
  int BETA = 1;

  // Get dvdlon = d(vecModFab[1])/dlon
  //            = d(vecModFab[1])/dalpha * dalpha/dlon +
  //              d(vecModFab[1])/beta * dbeta/dlon
  //            = diffVecModFab[0][1] * partialsFab[0] +
  //              diffVecModFab[1][1] * partialsFab[1]
  FArrayBox dvdlonFab(a_box, 1);
  workFab.copy(diffVecModFab[ALPHA], MERIDIONAL, 0);
  workFab.mult(partialsFab, ALPHALON, 0);
  dvdlonFab.copy(workFab);
  workFab.copy(diffVecModFab[BETA], MERIDIONAL, 0);
  workFab.mult(partialsFab, BETALON, 0);
  dvdlonFab += workFab;

  // Get dudlat = d(vecModFab[0])/dlat
  //            = d(vecModFab[0])/dalpha * dalpha/dlat +
  //              d(vecModFab[0])/beta * dbeta/dlat
  //            = diffVecModFab[0][0] * partialsFab[2] +
  //              diffVecModFab[1][0] * partialsFab[3]
  FArrayBox dudlatFab(a_box, 1);
  workFab.copy(diffVecModFab[ALPHA], ZONALCOSLAT, 0);
  workFab.mult(partialsFab, ALPHALAT, 0);
  dudlatFab.copy(workFab);
  workFab.copy(diffVecModFab[BETA], ZONALCOSLAT, 0);
  workFab.mult(partialsFab, BETALAT, 0);
  dudlatFab += workFab;

  // Set curlFab = (dvdlonFab - dudlatFab) / coslatFab.
  CH_assert(a_curlFab.box().contains(a_box));
  a_curlFab.copy(dvdlonFab);
  a_curlFab -= dudlatFab;
  a_curlFab.divide(coslatFab, 0, 0);
}

#include "NamespaceFooter.H"
