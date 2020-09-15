#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CubedSphereShellPanelCS.H"
#include "BoxIterator.H"
#include <cmath>
#include "MBMiscUtil.H"

#include "CubedSphereShellF_F.H"
#include "FourthOrderUtil.H"
#include "FourthOrderUtilF_F.H"
#include "FCDivergenceF_F.H"
#include "AdvectOpF_F.H"

#include "NamespaceHeader.H"
// #include "UsingNamespace.H"

CubedSphereShellPanelCS::CubedSphereShellPanelCS(
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
  m_height = 1;
  m_radius = 1;
}

CubedSphereShellPanelCS::~CubedSphereShellPanelCS()
{
}

void
CubedSphereShellPanelCS::setFlatMap(bool a_flatMap)
{
  m_flatMap = a_flatMap;
  if (m_flatMap) 
    m_realDim = SpaceDim;
  else
    m_realDim = 3;
}

RealVect
CubedSphereShellPanelCS::realCoord(const RealVect& a_Xi) const
{
  // Real coordinate cannot be stored in a RealVect
  MayDay::Error("Not implemented for CubedSphereShellPanelCS!!");
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
CubedSphereShellPanelCS::realCoord(FArrayBox& a_x, const FArrayBox& a_Xi,
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
#if CH_SPACEDIM >= 3
          a_x(iv, 2) = a_Xi(iv, 2) * m_height;
#endif
        }
    }
  else // !m_flatMap
    {
      if (m_verticalMap.isNull())
      {
        RealVect offset = centeringOffset(a_box, m_dx);
        FORT_CUBEDSPHERESHELLMAPPEDTOREAL(CHF_FRA(a_x),
                                        CHF_CONST_FRA(a_Xi),
                                        CHF_CONST_INT(m_nPanel),
                                        CHF_BOX(a_box),
                                        CHF_CONST_REALVECT(m_dx),
                                        CHF_CONST_REALVECT(offset),
                                        CHF_CONST_INTVECT(m_ix));
      }
      else
      { 
        // Pass in the 1D vertical mapping as an array
        const int npts = m_verticalMap->getNumPoints();
        const Real* pts = m_verticalMap->getPoints();
        CH_assert(pts);
        RealVect offset = centeringOffset(a_box, m_dx);
        // this handles horizontal
        FORT_CUBEDSPHERESHELLMAPPEDTOREAL(CHF_FRA(a_x),
                                        CHF_CONST_FRA(a_Xi),
                                        CHF_CONST_INT(m_nPanel),
                                        CHF_BOX(a_box),
                                        CHF_CONST_REALVECT(m_dx),
                                        CHF_CONST_REALVECT(offset),
                                        CHF_CONST_INTVECT(m_ix));
        // this handles vertical
        // TODO this is probably slow and horrible
        // but it passes the level exchange test
#if CH_SPACEDIM >= 3
      BoxIterator bit(a_box);
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          if ( (a_Xi(iv,2)>=0.) && (a_Xi(iv,2)<=1.) ) {
            a_x(iv, 2) = m_verticalMap->getValue(a_Xi(iv,2));
            a_x(iv, 2)*= m_height;
          } else { a_x(iv,2) = 0.;}
        }
#endif
        
    }
  }
}


// given coordinate in real space, return its location in the mapped space
RealVect
CubedSphereShellPanelCS::mappedCoord(const RealVect& a_x) const
{
  // Real coordinate cannot be stored in a RealVect
  pout() << "Not implemented for CubedSphereShellPanelCS!!" << endl;
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
CubedSphereShellPanelCS::dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const
{
  //pout() << "Not implemented for CubedSphereShellPanelCS!!" << endl;
  //CH_assert(false);

  CH_assert((m_nPanel >= 0) && (m_nPanel <= 5));

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
    // dx/dr
    else if (a_dirXi == 2)
    {
      switch (m_nPanel)
      {
        case 0:
          value = dInvDelta;
          break;

        case 1:
          value = - dGX * dInvDelta;
          break;

        case 2:
          value = - dInvDelta;
          break;

        case 3:
          value = dGX * dInvDelta;
          break;

        case 4:
          value = - dGY * dInvDelta;
          break;

        case 5:
          value = dGY * dInvDelta;
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
    // dy/dr
    else if (a_dirXi == 2)
    {
      switch (m_nPanel)
      {
        case 0:
          value = dGX * dInvDelta;
          break;

        case 1:
          value = dInvDelta;
          break;

        case 2:
          value = - dGX * dInvDelta;
          break;

        case 3:
          value = - dInvDelta;
          break;

        case 4:
          value = dGX * dInvDelta;
          break;

        case 5:
          value = dGX * dInvDelta;
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
    // dz/dr
    else if (a_dirXi == 2)
    {
      switch (m_nPanel)
      {
        case 0:
          value = dGY * dInvDelta;
          break;

        case 1:
          value = dGY * dInvDelta;
          break;

        case 2:
          value = dGY * dInvDelta;
          break;

        case 3:
          value = dGY * dInvDelta;
          break;

        case 4:
          value = dInvDelta;
          break;

        case 5:
          value = - dInvDelta;
          break;
      }
    }
  }

  if (a_dirXi != 2){
    Real r=m_radius;
#if CH_SPACEDIM>=3
    r+=a_Xi[2]*m_height;
#endif
    value*=r;}

  return value;
}

Real
CubedSphereShellPanelCS::dXidX(const RealVect& a_X, int a_dirX, int a_dirXi) const
{
  const Real& x=a_X[0];
  const Real& y=a_X[1];
  const Real& z=a_X[2];
  Real x2=x*x;
  Real y2=y*y;
  Real z2=z*z;

  Real value;

  if (a_dirX == 0){
    if (a_dirXi == 0){// dalpha/dx
      switch (m_nPanel){
      case 0: value = -y/(x2+y2); break;
      case 1: value = -y/(x2+y2); break;
      case 2: value = -y/(x2+y2); break;
      case 3: value = -y/(x2+y2); break;
      case 4: value = 0; break;
      case 5: value = 0; break;}}
    else if (a_dirXi == 1){ // dbeta/dx
      switch(m_nPanel){
      case 0: value = -z/(z2+x2); break;
      case 1: value = 0; break;
      case 2: value = z/(z2+x2); break;
      case 3: value = 0; break;
      case 4: value = -z/(z2+x2); break;
      case 5: value = -z/(z2+x2); break;}}
    else
      switch(m_nPanel){ // dr/dx
      case 0: value = x/sqrt(x2+y2+z2); break;
      case 1: value = x/sqrt(x2+y2+z2); break;
      case 2: value = x/sqrt(x2+y2+z2); break;
      case 3: value = x/sqrt(x2+y2+z2); break;
      case 4: value = x/sqrt(x2+y2+z2); break;
      case 5: value = x/sqrt(x2+y2+z2); break;}}
  else if (a_dirX == 1){
    if (a_dirXi == 0){// dalpha/dy
      switch (m_nPanel){
      case 0: value = x/(y2+x2); break;
      case 1: value = x/(y2+x2); break;
      case 2: value = x/(y2+x2); break;
      case 3: value = x/(y2+x2); break;
      case 4: value = z/(y2+z2); break;
      case 5: value = -z/(y2+z2); break;}}
    else if (a_dirXi == 1){ // dbeta/dy
      switch(m_nPanel){
      case 0: value = 0; break;
      case 1: value = -z/(y2+z2); break;
      case 2: value = 0; break;
      case 3: value = z/(y2+z2); break;
      case 4: value = 0; break;
      case 5: value = 0; break;}}
    else
      switch(m_nPanel){ // dr/dy
      case 0: value = y/sqrt(x2+y2+z2); break;
      case 1: value = y/sqrt(x2+y2+z2); break;
      case 2: value = y/sqrt(x2+y2+z2); break;
      case 3: value = y/sqrt(x2+y2+z2); break;
      case 4: value = y/sqrt(x2+y2+z2); break;
      case 5: value = y/sqrt(x2+y2+z2); break;}}
  else {
    if (a_dirXi == 0){// dalpha/dz
      switch (m_nPanel){
      case 0: value = 0; break;
      case 1: value = 0; break;
      case 2: value = 0; break;
      case 3: value = 0; break;
      case 4: value = -y/(y2+z2); break;
      case 5: value = y/(y2+z2); break;}}
    else if (a_dirXi == 1){ // dbeta/dz
      switch(m_nPanel){
      case 0: value = x/(z2+x2); break;
      case 1: value = y/(z2+y2); break;
      case 2: value = -x/(z2+x2); break;
      case 3: value = -y/(z2+y2); break;
      case 4: value = x/(z2+x2); break;
      case 5: value = x/(z2+x2); break;}}
    else
      switch(m_nPanel){ // dr/dz
      case 0: value = z/sqrt(x2+y2+z2); break;
      case 1: value = z/sqrt(x2+y2+z2); break;
      case 2: value = z/sqrt(x2+y2+z2); break;
      case 3: value = z/sqrt(x2+y2+z2); break;
      case 4: value = z/sqrt(x2+y2+z2); break;
      case 5: value = z/sqrt(x2+y2+z2); break;}}

  return value;
}

void
CubedSphereShellPanelCS::getNodeRealCoordinates(
                                             FArrayBox& a_nodeCoords,
                                             const Box& a_box) const
{
  if (m_verticalMap.isNull())
  {
    FORT_CUBEDSPHERESHELLNODEREALFORVIZ(CHF_FRA(a_nodeCoords),
                                        CHF_CONST_INT(m_nPanel),
                                        CHF_BOX(a_box),
                                        CHF_CONST_REALVECT(m_dx),
                                        CHF_CONST_INTVECT(m_ix));
  }
  else
  { 
    // Pass in the 1D vertical mapping as an array
    const int npts = m_verticalMap->getNumPoints();
    const Real* pts = m_verticalMap->getPoints();
    CH_assert(pts);
    FORT_CUBEDSPHERESHELLVERTMAPNODEREALFORVIZ(CHF_FRA(a_nodeCoords),
                                          CHF_CONST_INT(m_nPanel),
                                          CHF_BOX(a_box),
                                          CHF_CONST_R1D(pts,npts),
                                          CHF_CONST_REALVECT(m_dx),
                                          CHF_CONST_INTVECT(m_ix));
  }
}


void 
CubedSphereShellPanelCS::getFaceMappedCoordinates(
    FArrayBox& a_faceCoords,
    const int a_dir,
    const Box& a_box) const
{
  CH_assert(a_box.type() == IntVect::Unit - BASISV(a_dir));
  getMappedCoordinates(a_faceCoords, a_box);
}


void CubedSphereShellPanelCS::getCellMappedCoordinates(
                                                    FArrayBox& a_cellCoords,
                                                    const Box& a_box
                                                    ) const
{
  CH_assert(a_box.ixType() == IndexType::TheCellType());
  getMappedCoordinates(a_cellCoords, a_box);
}


void CubedSphereShellPanelCS::getNodeMappedCoordinates(
                                                    FArrayBox& a_nodeCoords,
                                                    const Box& a_box
                                                    ) const
{
  CH_assert(a_box.ixType() == IndexType::TheNodeType());
  getMappedCoordinates(a_nodeCoords, a_box);
}


void CubedSphereShellPanelCS::getMappedCoordinates(
                                                   FArrayBox& a_coords,
                                                   const Box& a_box
                                                   ) const
{
  CH_assert(a_coords.box().contains(a_box));
  RealVect offset = centeringOffset(a_box, m_dx);
  FORT_CUBEDSPHERESHELLMAPPEDCOORDS(CHF_FRA(a_coords),
                                    CHF_CONST_INT(m_nPanel),
                                    CHF_BOX(a_box),
                                    CHF_CONST_REALVECT(m_dx),
                                    CHF_CONST_REALVECT(offset),
                                    CHF_CONST_INTVECT(m_ix));
}


//-----------------------------------------------------------------------
void
CubedSphereShellPanelCS::
volFlux(FluxBox& a_volFlux,
        const FluxBox& a_Nt,
        const Box& a_box) const
{
  CH_TIME("CubedSphereShellPanelCS::volFlux");

  // Note that X needs to have one more ghost cell
  Box bx1 = grow(a_box, 1);
  CH_assert(a_Nt.box().contains(bx1));
  FluxBox NtX(bx1, SpaceDim);

  int compVol0 = m_volInterval.begin();
  int compVol1 = m_volInterval.end();

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
    FArrayBox& thisNtXdir = NtX[dir];

    const Box& thisNtXBox = thisNtXdir.box();
    RealVect offset = centeringOffset(thisNtXBox, m_dx);

    BoxIterator bit(thisNtXBox);
    // this is going to be slow, but we can
    // eventually move this into fortran
    for (bit.begin(); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        RealVect mappedLoc = m_dx*iv + offset;

        // int signX = (sin(mappedLoc[0] + 0.25 * M_PI) < 0.) ? -1 : 1;
        Real dGX = signX*tan(mappedLoc[0] - 0.25 * M_PI);
        // int signY = ((mappedLoc[0] + 0.25 * M_PI) > (2 * M_PI)) ? -1 : 1;
        Real dGY = signY*tan(mappedLoc[1] - 0.25 * M_PI);
        Real dDelta = sqrt(1.0 + dGX*dGX + dGY*dGY);
        thisNtXdir(iv, compVol0) = 0.5 * dGX / dDelta;
        thisNtXdir(iv, compVol1) = 0.5 * dGY / dDelta;
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
CubedSphereShellPanelCS::getLineMetrics(FArrayBox& a_metrics,
                              const Box& a_box) const
{
  CH_TIME("CubedSphereShellPanelCS::getLineMetrics");

  CH_assert(a_metrics.nComp() == SpaceDim);
  RealVect offset = centeringOffset(a_box, m_dx);
  // NOTE: height and radius are ignored in this routine!
  FORT_CUBEDSPHERESHELLLINEMETRICS(CHF_FRA(a_metrics),
                                   CHF_CONST_INT(m_nPanel),
                                   CHF_BOX(a_box),
                                   CHF_CONST_REALVECT(m_dx),
                                   CHF_CONST_REAL(m_height),
                                   CHF_CONST_REAL(m_radius),
                                   CHF_CONST_REALVECT(offset),
                                   CHF_CONST_INTVECT(m_ix));
}


void
CubedSphereShellPanelCS::cellVol(FArrayBox& a_vol,
                              const FluxBox& a_N,
                              const Box& a_box) const
{
  CH_TIME("CubedSphereShellPanelCS::cellVol");

  // Pass in the 1D vertical mapping as an array
#if CH_SPACEDIM>=3
  CH_assert(!m_verticalMap.isNull());
  const int npts = m_verticalMap->getNumPoints();
  const Real* pts = m_verticalMap->getPoints();
  CH_assert(pts);
#else
  const int npts = 0;
  const Real* pts = NULL;
#endif

  FORT_CUBEDSPHERESHELLCELLVOL(CHF_FRA1(a_vol, 0),
                               CHF_CONST_INT(m_nPanel),
                               CHF_BOX(a_box),
                               CHF_CONST_R1D(pts,npts),
                               CHF_CONST_REALVECT(m_dx),
                               CHF_CONST_REAL(m_height),
                               CHF_CONST_REAL(m_radius),
                               CHF_CONST_INTVECT(m_ix));
}

void CubedSphereShellPanelCS::getN(FluxBox& a_N, const Box& a_box) const
{
  CH_TIME("CubedSphereShellPanelCS::getN");

  a_N.setVal(0.0);
  // // a_N[0] and a_N[1] both have 2*2 = 4 components

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      FArrayBox& Nfab = a_N[idir];
      // added by petermc, 12 Apr 2011
      Box faceBox(a_box);
      faceBox.surroundingNodes(idir);
      // In 2D:
      // If idir == 0, area of face is length of edge, which is dx[0].
      // If idir == 1, area of face is length of edge, which is dx[1].
      // In 3D, each of these face areas is multiplied by dx[2],
      // and if idir == 2, which is new, area of face is always dx[0]*dx[1].
      FORT_CUBEDSPHERESHELLNORMAL(CHF_FRA(Nfab),
                                  CHF_CONST_INT(m_nPanel),
                                  CHF_BOX(faceBox),
                                  CHF_CONST_INT(idir),
                                  CHF_CONST_REALVECT(m_dx),
                                  CHF_CONST_INTVECT(m_ix));
      // Kludge added by petermc, 27 Sep 2012.
      // Re-kludged to add height scaling by hans, 20 Oct 2012.
      if (SpaceDim == 3){
        if (idir < 2){
          // height of face in physical space
          Real dhReal = m_height * m_dx[2];
          // In 2D, the components of Nfab are:
          // 0: dx[0]/dxi[0] the ONLY nonzero when idir == 0
          // 1: dx[0]/dxi[1]
          // 2: dx[1]/dxi[0]
          // 3: dx[1]/dxi[1] the ONLY nonzero when idir == 1
          // In 3D, the components of Nfab are:
          // 0: dx[0]/dxi[0] == old[0]*dhReal
          // 1: dx[0]/dxi[1] == old[1]*dhReal
          // 2: dx[0]/dxi[2] == 0.
          // 3: dx[1]/dxi[0] == old[2]*dhReal
          // 4: dx[1]/dxi[1] == old[3]*dhReal
          // 5: dx[1]/dxi[2] == 0.
          // 6: dx[2]/dxi[0] == 0.
          // 7: dx[2]/dxi[1] == 0.
          // 8: dx[2]/dxi[2] == 0.
          // So we must set:
          // new[4]=old[3]
          // new[3]=old[2];
          // new[2]=0.;
          // and then multiply new[0], new[1], new[3], new[4] by dhReal.
          Nfab.copy(Nfab, 3, 4); // Nfab[4] := Nfab[3]
          Nfab.copy(Nfab, 2, 3); // Nfab[3] := Nfab[2]
          Nfab.setVal(0., 2); // Nfab[2] := 0.
          Nfab.mult(dhReal, 0); // Nfab[0] *= dhReal
          Nfab.mult(dhReal, 1); // Nfab[1] *= dhReal
          Nfab.mult(dhReal, 3); // Nfab[3] *= dhReal
          Nfab.mult(dhReal, 4); // Nfab[4] *= dhReal
        }else{
          Nfab.mult(m_height,8);
        }
      }
    }

  a_N*=m_radius*m_radius;
}

/// computes cell-averaged J
void CubedSphereShellPanelCS::getAvgJ(FArrayBox& a_avgJ,
                                   const Box& a_box) const
{
  CH_TIME("CubedSphereShellPanelCS::getAvgJ");

  // we don't need N the way we're computing this, but it's required by
  // the API
  FluxBox Nbogus;

  cellVol(a_avgJ, Nbogus, a_box);

  Real mappedVol = m_dx.product();
  a_avgJ /= mappedVol;
}

/// computes cell-averaged J
void CubedSphereShellPanelCS::getAvgJ(FArrayBox& a_avgJ,
                                   const FluxBox& a_volFlux,
                                   const Box& a_box) const
{
  MayDay::Error("This shouldn't be used!");

  a_avgJ.setVal(0.0, a_box, 0, a_avgJ.nComp());
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

  for (int dir=0; dir != SpaceDim; ++dir)
    {
      FORT_FCDIVERGENCE(CHF_CONST_FRA(a_volFlux[dir]),
                        CHF_FRA(a_avgJ),
                        CHF_BOX(a_box),
                        CHF_CONST_REAL(panelDxSign[dir]),
                        CHF_INT(dir));
    }
}

/// computes cell-averaged 1/J
void CubedSphereShellPanelCS::getAvgJinverse(
                                          FluxBox& a_avgJinverse,
                                          const Box& a_box
                                          ) const
{
  pout() << "Not implemented for CubedSphereShellPanelCS!!" << endl;
  CH_assert(false);


  for (int idir = 0; idir < SpaceDim; idir++)
    {
      FArrayBox& avgJinvFab = a_avgJinverse[idir];
      RealVect offset = centeringOffset(avgJinvFab.box(), m_dx);
      CH_assert(SpaceDim == 2);
      FORT_CUBEDSPHERESHELLAVGJINV(CHF_FRA1(avgJinvFab, 0),
                                   CHF_CONST_INT(m_nPanel),
                                   CHF_BOX(a_box),
                                   CHF_CONST_INT(idir),
                                   CHF_CONST_REALVECT(m_dx),
                                   CHF_CONST_REALVECT(offset),
                                   CHF_CONST_INTVECT(m_ix));
    }
}

/// Jacobian evaluated at location X in mapped space
Real CubedSphereShellPanelCS::pointwiseJ(const RealVect& a_Xi) const
{
  pout() << "Not implemented for CubedSphereShellPanelCS!!" << endl;
  CH_assert(false);


  Real dGX = tan(a_Xi[0] - 0.25 * M_PI);
  Real dGY = tan(a_Xi[1] - 0.25 * M_PI);

  Real dDelta = sqrt(1.0 + dGX * dGX + dGY * dGY);

  return (1.0 + dGX * dGX) * (1.0 + dGY * dGY) / (dDelta * dDelta * dDelta);
}

/// Jacobian evaluated at locations Xi in mapped space
void CubedSphereShellPanelCS::pointwiseJ(
                                      FArrayBox& a_J,
                                      const FArrayBox& a_Xi,
                                      const Box& a_box
                                      ) const
{
  CH_TIME("CubedSphereShellPanelCS::pointwiseJ");

  FORT_CUBEDSPHERESHELLPOINTJ(CHF_FRA1(a_J, 0),
                              CHF_CONST_FRA(a_Xi),
                              CHF_CONST_INT(m_nPanel),
                              CHF_BOX(a_box),
                              CHF_CONST_REALVECT(m_dx),
                              CHF_CONST_INTVECT(m_ix),
                              CHF_CONST_REAL(m_height),
                              CHF_CONST_REAL(m_radius));
}

/// returns integral of divergence over mapped-grid cells
void CubedSphereShellPanelCS::computeDivergence(
                                             FArrayBox& a_divF,
                                             const FluxBox& a_F,
                                             const FluxBox& a_N,
                                             const Box& a_box,
                                             Interval& divInterval
) const
{
  pout() << "Not implemented for CubedSphereShellPanelCS!!" << endl;
  CH_assert(false);

}

/// transform a point from mapped-coordinate basis to latitude-longitude
/// coordinate basis
void CubedSphereShellPanelCS::pointTransformEquiangularToLonLat(
  const RealVect& a_xi,
  RealVect& a_rllXi
) const
{
  Real eps = DBL_EPSILON;
  FORT_CUBEDSPHERESHELLEQUIANGULARTOLONLAT(CHF_REALVECT(a_rllXi),
                                        CHF_CONST_REALVECT(a_xi),
                                        CHF_CONST_INT(m_nPanel),
                                        CHF_CONST_REAL(eps));

  // Radial component is mapped and scaled by height
  if (!m_verticalMap.isNull()){
    Real xi = a_xi[2];
    if (xi<0.) xi=0.; // these are clunky and only here because levels
    if (xi>1.) xi=1.; // make ghost cells in vertical for no clear reason
    a_rllXi[2] = m_verticalMap->getValue(xi);
    a_rllXi[2] *= m_height;
  }
}


/// transform a point from mapped-coordinate basis to Cartesian
/// coordinate basis
void
CubedSphereShellPanelCS::pointTransformEquiangularToCartesian(
                                                           const RealVect& a_xi,
                                                           Real * a_xyz)
  const
{
  //  pout() << "Not implemented for CubedSphereShellPanelCS!!" << endl;
  //  CH_assert(false);

  //  CH_assert(SpaceDim == 2);
  FORT_CUBEDSPHERESHELLEQUIANGULARTOCARTESIAN(CHF_R1D(a_xyz, 3),
                                              CHF_CONST_REALVECT(a_xi),
                                              CHF_CONST_INT(m_nPanel),
                                              CHF_CONST_REAL(m_radius),
                                              CHF_CONST_REAL(m_height));
}


/// transform a FAB of points from mapped-coordinate basis to latitude-longitude
/// coordinate basis
void CubedSphereShellPanelCS::fabTransformEquiangularToLonLat(
                                                           const FArrayBox& a_xiFab,
                                                           FArrayBox& a_rllXiFab
                                                           ) const
{
  CH_TIME("CubedSphereShellPanelCS::fabTransformEquiangularToLonLat");

  CH_assert(a_xiFab.box() == a_rllXiFab.box());
  Real eps = DBL_EPSILON;
  // convert horizontal coordinates
  FORT_CUBEDSPHERESHELLFABEQUIANGULARTOLONLAT(CHF_FRA(a_rllXiFab),
                                              CHF_CONST_FRA(a_xiFab),
                                              CHF_CONST_INT(m_nPanel),
                                              CHF_CONST_REAL(m_height),
                                              CHF_CONST_REAL(eps));
#if CH_SPACEDIM >= 3
  // convert vertical coordinates
  // TODO this is probably a slow mess
  if (!m_verticalMap.isNull()){
  Box box = a_xiFab.box();
  BoxIterator bit(box);
  for (bit.begin();bit.ok();++bit){
    IntVect iv = bit();
    Real xi = a_xiFab(iv,2);
    if (xi<0.) xi=0.; // these are clunky and only here because levels
    if (xi>1.) xi=1.; // make ghost cells in vertical for no clear reason
    a_rllXiFab(iv,2)  = m_verticalMap->getValue(xi);
    a_rllXiFab(iv,2) *= m_height;
  }
  }
#endif
}


/// transform a FAB of points from mapped-coordinate basis to Cartesian
/// coordinate basis
void
CubedSphereShellPanelCS::fabTransformEquiangularToCartesian(
                                                         const FArrayBox& a_xiFab,
                                                         FArrayBox& a_xyzFab) const
{
  CH_TIME("CubedSphereShellPanelCS::fabTransformEquiangularToCartesian");

  CH_assert(a_xiFab.box() == a_xyzFab.box());
  FORT_CUBEDSPHERESHELLFABEQUIANGULARTOCARTESIAN(CHF_FRA(a_xyzFab),
                                                 CHF_CONST_FRA(a_xiFab),
                                                 CHF_CONST_INT(m_nPanel),
                                                 CHF_CONST_REAL(m_radius),
                                                 CHF_CONST_REAL(m_height));
}


/// transform a vector from mapped-coordinate basis to real-coordinate basis
void CubedSphereShellPanelCS::vectorTransformEquiangularToCartesian(
  const RealVect& a_xi,
  const Real * a_vecCS,
  Real * a_vecXYZ
) const
{
    for (int i=0;i<3;i++){
      a_vecXYZ[i]=0;
      for (int j=0;j<SpaceDim;j++)
        a_vecXYZ[i]+=a_vecCS[j]*dXdXi(a_xi,i,j);}
}

/// transform a vector from real-coordinate basis to mapped-coordinate basis
void CubedSphereShellPanelCS::vectorTransformCartesianToEquiangular(
  const RealVect& a_xi,
  const Real * a_vecXYZ,
  Real * a_vecCS
) const
{
  RealVect x;
  pointTransformEquiangularToCartesian(a_xi,(Real*)&x);

  for (int i=0;i<SpaceDim;i++){
    a_vecCS[i]=0;
    for (int j=0;j<3;j++)
      a_vecCS[i]+=a_vecXYZ[j]*dXidX(x,j,i);}
}

/// transform a vector from latitude-longitude coordinate basis to
/// mapped-coordinate basis
void CubedSphereShellPanelCS::vectorTransformLatLonToEquiangular(
  const RealVect& a_xi,
  const RealVect& a_vecRLL,
  RealVect& a_vecCS
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

#if CH_SPACEDIM >= 3
  // vertical
  if (m_verticalMap.isNull())
    {
      a_vecCS[2] = a_vecRLL[2];
    }
  else
    {
      Real dzdxi = m_verticalMap->getDerivative(a_xi[2]);
      a_vecCS[2] = a_vecRLL[2] / (dzdxi*m_height);
    }
#endif

  // horizontal (lat lon)
  switch(m_nPanel)
  {
    // Equatorial panels
    case 0:
    case 1:
    case 2:
    case 3:
      lat = atan(dGY / sqrt(1.0 + dGX * dGX));
      vecLon = a_vecRLL[0] / cos(lat); // why do we do this?

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
      if (gnoRadius==0.){ // on the pole
        // TODO is this right?
        // either way, need to propogate to the FAB version of this
        a_vecCS[0] = vecLon;
        a_vecCS[1] = vecLon;
        break;
      } 

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
       if (gnoRadius==0.){ // on the pole
         // TODO is this right?
         a_vecCS[0] = vecLon;
         a_vecCS[1] = vecLon;
         break;
       } 

      a_vecCS[0] =
        dGY / (1.0 + dGX * dGX) * vecLon
        + dDelta2 * dGX / ((1.0 + dGX * dGX) * gnoRadius) * vecLat;

      a_vecCS[1] =
        - dGX / (1.0 + dGY * dGY) * vecLon
        + dDelta2 * dGY / ((1.0 + dGY * dGY) * gnoRadius) * vecLat;
      break;
  }

  // TODO copied from the fortran routine for Fabs, why do we do this?
  a_vecCS[0] /= m_radius;
  a_vecCS[1] /= m_radius;
}

/// transform a FAB of vectors from latitude-longitude coordinate basis to
/// mapped-coordinate basis
void CubedSphereShellPanelCS::fabVectorTransformLatLonToEquiangular(
    const FArrayBox& a_xiFab,
    const FArrayBox& a_vecRLLFab,
    FArrayBox& a_vecCSFab) const
{
  CH_assert((m_nPanel >= 0) && (m_nPanel <= 5));
  //CH_assert(a_xiFab.box() == a_vecRLLFab.box());
  //CH_assert(a_vecCSFab.box().contains(a_xiFab.box()));
  // this is horizontal only
  FORT_CUBEDSPHERESHELLFABVECTORLATLONTOEQUIANGULAR(
      CHF_CONST_FRA(a_xiFab),
      CHF_CONST_FRA(a_vecRLLFab),
      CHF_FRA(a_vecCSFab),
      CHF_CONST_REAL(m_height),
      CHF_CONST_REAL(m_radius),
      CHF_CONST_INT(m_nPanel),
      CHF_BOX(a_vecCSFab.box()));

#if CH_SPACEDIM >= 3
  // convert in vertical as well
  // TODO this is probably really slow
  if (!m_verticalMap.isNull()){
  BoxIterator bit(a_vecCSFab.box());
  for(bit.begin();bit.ok();++bit){
    IntVect iv = bit();

    Real dzdxi = m_verticalMap->getDerivative(a_xiFab(iv,2));
    a_vecCSFab(iv,2) = a_vecRLLFab(iv,2) / (m_height*dzdxi);
  }
  }
#endif
}


/// transform a FAB of vectors from mapped-coordinate basis to
/// latitude-longitude coordinate basis
void CubedSphereShellPanelCS::fabVectorTransformEquiangularToLatLon(
                                                              const FArrayBox& a_xiFab,
                                                              const FArrayBox& a_vecCSFab,
                                                              FArrayBox& a_vecRLLFab) const
{
  CH_assert((m_nPanel >= 0) && (m_nPanel <= 5));
  //CH_assert(a_xiFab.box() == a_vecRLLFab.box());
  //CH_assert(a_xiFab.box() == a_vecCSFab.box());
  //CH_assert(SpaceDim == 2);
  Box box=a_vecRLLFab.box();
  FORT_CUBEDSPHERESHELLFABVECTOREQUIANGULARTOLATLON(CHF_CONST_FRA(a_xiFab),
                                                    CHF_CONST_FRA(a_vecCSFab),
                                                    CHF_FRA(a_vecRLLFab),
                                                    CHF_CONST_REAL(m_height),
                                                    CHF_CONST_REAL(m_radius),
                                                    CHF_CONST_INT(m_nPanel),
                                                    CHF_BOX(box));
}


/// transform a FAB of SpaceDim-vectors from mapped-coordinate basis to
/// real-coordinate basis at cell centers
void CubedSphereShellPanelCS::vectorTransformMappedToRealCenterFab(
                                          FArrayBox& a_vectorFab
) const
{
  Real vecCS[3];
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

      vecCS[0] = a_vectorFab(iv,0);
      vecCS[1] = a_vectorFab(iv,1);
#if CH_SPACEDIM>=3
      vecCS[2] = a_vectorFab(iv,2);
#endif

      vectorTransformEquiangularToCartesian(ABcenter, vecCS, vecXYZ);

      a_vectorFab(iv,0) = vecXYZ[0];
      a_vectorFab(iv,1) = vecXYZ[1];
      a_vectorFab(iv,2) = vecXYZ[2];
    }
}

/// transform a FAB of SpaceDim-vectors from real-coordinate basis to
/// mapped-coordinate basis at cell centers
void
CubedSphereShellPanelCS::vectorTransformRealToMappedCenterFab(
                                                           FArrayBox& a_vectorFab
) const
{
  CH_assert(a_vectorFab.nComp() == 3); // use 3 comps on input, SpaceDim on output
  const Box& bx = a_vectorFab.box();
  FORT_CUBEDSPHERESHELLVECTORREALTOMAPPEDCENTER(CHF_FRA(a_vectorFab),
                                             CHF_CONST_INT(m_nPanel),
                                             CHF_BOX(bx),
                                             CHF_CONST_REALVECT(m_dx),
                                             CHF_CONST_INTVECT(m_ix));
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
CubedSphereShellPanelCS::
getNMatrixEntry(const RealVect& a_Xi,
                int a_s, int a_d, int a_d1,
                int a_row, int a_column) const
{
  pout() << "Not implemented for CubedSphereShellPanelCS!!" << endl;
  CH_assert(false);


  Real entry = 0.0;
  // simplest thing -- if row = s, then A_ij = delta_dj
  if (a_row == a_s)
  {
    if (a_column == a_d) entry = 1.0;
  }
  else if (a_column == a_d1)
  {
    // We can't use realCoord() here, since that method assumes that the
    // dimensions of the domain and the range spaces are the same.
    // Instead, we directly compute the (a_row)th component of X.
    RealVect ll;
    pointTransformEquiangularToLonLat(a_Xi, ll);
    entry = ll[a_row];
  }
  else
  {
    entry = dXdXi(a_Xi, a_row, a_column);
  }

  return entry;
}

//-----------------------------------------------------------------------
Real
CubedSphereShellPanelCS::getN(const RealVect& a_Xi,
                           int a_s, int a_d, int a_d1) const
{
  pout() << "Not implemented for CubedSphereShellPanelCS!!" << endl;
  CH_assert(false);


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
          MayDay::Error("CubedSphereShellPanelCS::getN dimension must be 0 or 1");
        }
    }

  return retval;
}

//-----------------------------------------------------------------------

void
CubedSphereShellPanelCS::contravariantMetric(FArrayBox& a_metric,
                                          int a_dir) const
{
  CH_TIME("CubedSphereShellPanelCS::contravariantMetric");

  //  pout() << "Not implemented for CubedSphereShellPanelCS!!" << endl;
  //  CH_assert(false);

  //  MayDay::Warning("Not sure CubedSphereShellPanelCS::contravariantMetric is correct!!");

  CH_assert(a_metric.nComp() == SpaceDim);
  const Box& bx = a_metric.box();

  FArrayBox xiFab(bx, SpaceDim);
  getCenterMappedCoordinates(xiFab, bx); // according to centering of bx

  FORT_CUBEDSPHERESHELLCONTRAVARIANTMETRIC(CHF_CONST_FRA(xiFab),
                                           CHF_FRA(a_metric),
                                           CHF_CONST_INT(a_dir),
                                           CHF_CONST_REAL(m_height),
                                           CHF_CONST_REAL(m_radius));
}

//-----------------------------------------------------------------------
void
CubedSphereShellPanelCS::orthonormalize(
                                     const FArrayBox& a_csFab,
                                     FArrayBox& a_orthoFab,
                                     const Box& a_box,
                                     int a_idir,
                                     const IntVect& a_csComps,
                                     const IntVect& a_orthoComps) const
{
  CH_TIME("CubedSphereShellPanelCS::orthonormalize");

  FArrayBox xiFab(a_box, SpaceDim);
  getCenterMappedCoordinates(xiFab, a_box); // according to centering of a_box

  FORT_CUBEDSPHERESHELLORTHONORMALIZE(CHF_CONST_FRA(xiFab),
                                      CHF_CONST_FRA(a_csFab),
                                      CHF_FRA(a_orthoFab),
                                      CHF_CONST_INT(a_idir),
                                      CHF_CONST_INTVECT(a_csComps),
                                      CHF_CONST_INTVECT(a_orthoComps));
}

//-----------------------------------------------------------------------
void
CubedSphereShellPanelCS::deorthonormalize(
                                       const FArrayBox& a_orthoFab,
                                       FArrayBox& a_csFab,
                                       const Box& a_box,
                                       int a_idir,
                                       const IntVect& a_orthoComps,
                                       const IntVect& a_csComps) const
{
  CH_TIME("CubedSphereShellPanelCS::deorthonormalize");

  FArrayBox xiFab(a_box, SpaceDim);
  getCenterMappedCoordinates(xiFab, a_box); // according to centering of a_box

  FORT_CUBEDSPHERESHELLDEORTHONORMALIZE(CHF_CONST_FRA(xiFab),
                                        CHF_CONST_FRA(a_orthoFab),
                                        CHF_FRA(a_csFab),
                                        CHF_CONST_INT(a_idir),
                                        CHF_CONST_INTVECT(a_orthoComps),
                                        CHF_CONST_INTVECT(a_csComps));
}


//-----------------------------------------------------------------------
void
CubedSphereShellPanelCS::orthonormalizeVectorFluxes(
                                                    FArrayBox& a_fluxFab,
                                                    const Box& a_box,
                                                    const Interval& a_vectorIntv) const
{
  CH_TIME("CubedSphereShellPanelCS::orthonormalizeVectorFluxes");

  CH_assert(a_vectorIntv.size() == SpaceDim);
  IntVect baseComps = IntVect(D_DECL6(0, 1, 2, 3, 4, 5));
  FArrayBox vectorFluxFab(a_vectorIntv, a_fluxFab); // alias

  int idir = faceDimension(a_box);

  FArrayBox xiFab(a_box, SpaceDim);
  getCenterMappedCoordinates(xiFab, a_box); // according to centering of a_box

  FArrayBox orthoFab(a_box, SpaceDim);
  FORT_CUBEDSPHERESHELLORTHONORMALIZE(CHF_CONST_FRA(xiFab),
                                      CHF_CONST_FRA(vectorFluxFab),
                                      CHF_FRA(orthoFab),
                                      CHF_CONST_INT(idir),
                                      CHF_CONST_INTVECT(baseComps),
                                      CHF_CONST_INTVECT(baseComps));
  vectorFluxFab.copy(orthoFab);
  FORT_CUBEDSPHERESHELLORTHONORMALIZEFRAME(CHF_CONST_FRA(xiFab),
                                           CHF_FRA(vectorFluxFab),
                                           CHF_CONST_INT(idir));
}


//-----------------------------------------------------------------------
void
CubedSphereShellPanelCS::deorthonormalizeVectorFluxes(
                                                      FArrayBox& a_fluxFab,
                                                      const Box& a_box,
                                                      const Interval& a_vectorIntv) const
{
  CH_TIME("CubedSphereShellPanelCS::deorthornomalizeVectorFluxes");

  CH_assert(a_vectorIntv.size() == SpaceDim);
  IntVect baseComps = IntVect(D_DECL6(0, 1, 2, 3, 4, 5));
  FArrayBox vectorFluxFab(a_vectorIntv, a_fluxFab); // alias

  int idir = faceDimension(a_box);

  FArrayBox xiFab(a_box, SpaceDim);
  getCenterMappedCoordinates(xiFab, a_box); // according to centering of a_box

  FArrayBox deorthoFab(a_box, SpaceDim);
  FORT_CUBEDSPHERESHELLDEORTHONORMALIZE(CHF_CONST_FRA(xiFab),
                                        CHF_CONST_FRA(vectorFluxFab),
                                        CHF_FRA(deorthoFab),
                                        CHF_CONST_INT(idir),
                                        CHF_CONST_INTVECT(baseComps),
                                        CHF_CONST_INTVECT(baseComps));
  vectorFluxFab.copy(deorthoFab);
  FORT_CUBEDSPHERESHELLDEORTHONORMALIZEFRAME(CHF_CONST_FRA(xiFab),
                                             CHF_FRA(vectorFluxFab),
                                             CHF_CONST_INT(idir));
}


//-----------------------------------------------------------------------
void
CubedSphereShellPanelCS::orthonormalizeVectorFluxes(
                                                    FluxBox& a_flux,
                                                    const Interval& a_vectorIntv) const
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      FArrayBox& fluxFab = a_flux[idir];
      orthonormalizeVectorFluxes(fluxFab, fluxFab.box(), a_vectorIntv);
    }
}


//-----------------------------------------------------------------------
void
CubedSphereShellPanelCS::deorthonormalizeVectorFluxes(
                                                      FluxBox& a_flux,
                                                      const Interval& a_vectorIntv) const
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      FArrayBox& fluxFab = a_flux[idir];
      deorthonormalizeVectorFluxes(fluxFab, fluxFab.box(), a_vectorIntv);
    }
}


#if 0
//-----------------------------------------------------------------------

void
CubedSphereShellPanelCS::getOrthonormalizingMatrix(
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
CubedSphereShellPanelCS::getDeorthonormalizingMatrix(
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

#endif

//-----------------------------------------------------------------------

void
CubedSphereShellPanelCS::curlRadial(
                                    const FArrayBox& a_vecFab,
                                    FArrayBox& a_curlFab,
                                    const Box& a_box) const
{
  Box bx2 = grow(a_box, 2);
  CH_assert(a_vecFab.box().contains(bx2));
  FArrayBox xiFab(bx2, SpaceDim);
  getCenterMappedCoordinates(xiFab, bx2); // according to centering of a_box

  FArrayBox vecCenUnitFab(bx2, SpaceDim);
  // Convert vector components from natural basis to unit basis.
  FORT_CUBEDSPHERESHELLVECTORTOUNITBASIS(CHF_CONST_FRA(xiFab),
                                         CHF_CONST_FRA(a_vecFab),
                                         CHF_FRA(vecCenUnitFab),
                                         CHF_BOX(bx2));

  Tuple<FArrayBox, SpaceDim> gradVecFab;
  for (int idir=0; idir<SpaceDim; idir++)
    {
      gradVecFab[idir].define(a_box, SpaceDim);
      FORT_CENTEREDGRADIENT4THORDER(CHF_FRA(gradVecFab[idir]),
                                    CHF_CONST_FRA1(vecCenUnitFab, idir),
                                    CHF_CONST_REALVECT(m_dx),
                                    CHF_BOX(a_box));
    }
  CH_assert(a_curlFab.box().contains(a_box));
  FORT_CUBEDSPHERESHELLCURLR(CHF_CONST_FRA(xiFab),
                             CHF_CONST_FRA(gradVecFab[0]),
                             CHF_CONST_FRA(gradVecFab[1]),
                             CHF_FRA1(a_curlFab, 0),
                             CHF_BOX(a_box));
}

//-----------------------------------------------------------------------

void
CubedSphereShellPanelCS::curlSpherical(
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
  FORT_CUBEDSPHERESHELLPARTIALS(CHF_CONST_FRA(lonlatFab),
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

//-----------------------------------------------------------------------

void
CubedSphereShellPanelCS::divSpherical(
                                      const FArrayBox& a_vecFab,
                                      FluxBox& a_divFlux,
                                      const Box& a_box) const
{
  Box bx1 = grow(a_box, 1); // CELL-centered
  CH_assert(a_vecFab.box().contains(bx1));
  CH_assert(a_divFlux.box().contains(a_box));
  FArrayBox xiFab(bx1, SpaceDim);
  getCenterMappedCoordinates(xiFab, bx1);

  FArrayBox deltaFab(bx1, 1);
  FORT_CUBEDSPHERESHELLDELTA(CHF_CONST_FRA(xiFab),
                             CHF_FRA1(deltaFab, 0));
  FArrayBox vecModFab(bx1, SpaceDim);
  int ALPHA = 0;
  int BETA = 1;
  vecModFab.copy(a_vecFab);
  vecModFab.divide(deltaFab, 0, ALPHA);
  vecModFab.divide(deltaFab, 0, BETA);
#if (CH_SPACEDIM > 2)
  // vecModFab[RADIAL] = r^2 * vecFab[RADIAL]
  int RADIAL = 2;
  FArrayBox rFab(bx1, 1);
  rFab.copy(xiFab, RADIAL, 0);
  rFab += 1.; // radius is (0:1) + 1.
  vecModFab.mult(rFab, 0, RADIAL);
  vecModFab.mult(rFab, 0, RADIAL);
#endif

  for (int idirFace=0; idirFace<SpaceDim; idirFace++)
    {
      Box bxFace(a_box);
      bxFace.surroundingNodes(idirFace);

      // gradVecModFab will contain
      // gradient of vecModFab[0], then
      // gradient of vecModFab[1], then
      // gradient of vecModFab[2].
      FArrayBox gradVecModFab(bxFace, SpaceDim*SpaceDim);
      FORT_SECONDORDERGRADIENT(CHF_FRA(gradVecModFab),
                               CHF_CONST_FRA(vecModFab),
                               CHF_BOX(bxFace),
                               CHF_CONST_REALVECT(m_dx),
                               CHF_CONST_INT(idirFace));

      FArrayBox xiFaceFab(bxFace, SpaceDim);
      getCenterMappedCoordinates(xiFaceFab, bxFace);

      FArrayBox& divFab = a_divFlux[idirFace];
      FORT_CUBEDSPHERESHELLDIVONFACES(CHF_FRA1(divFab, 0),
                                      CHF_CONST_FRA(gradVecModFab),
                                      CHF_CONST_FRA(xiFaceFab),
                                      CHF_BOX(bxFace),
                                      CHF_CONST_INT(idirFace));
    }
}

//-----------------------------------------------------------------------

void
CubedSphereShellPanelCS::gradient(
                                  const FArrayBox& a_fab,
                                  FArrayBox& a_gradFab,
                                  const Box& a_box) const
{
  Box bx1 = grow(a_box, 1); // CELL-centered
  CH_assert(a_fab.box().contains(bx1));
  CH_assert(a_gradFab.box().contains(a_box));
  FArrayBox xiFab(bx1, SpaceDim);
  getCenterMappedCoordinates(xiFab, bx1);

  int ncomp = a_fab.nComp();
  // Order of components in diffFab and a_gradFab will be:
  // ncomp differences with dimension 0,
  // ncomp differences with dimension 1, etc.
  FArrayBox diffFab(a_box, SpaceDim*ncomp);
  for (int idir=0; idir<SpaceDim; idir++)
    {
      int baseComp = idir*ncomp;
      Interval intvl(baseComp, baseComp + ncomp-1);
      FArrayBox diffDirFab(intvl, diffFab);
      // FORT_UDIVCENTERGRAD returns undivided centered differences.
      FORT_UDIVCENTERGRAD(CHF_FRA(diffDirFab),
                          CHF_CONST_FRA(a_fab),
                          CHF_BOX(a_box),
                          CHF_CONST_INT(idir));
      diffDirFab /= m_dx[idir];
    }
  FORT_CUBEDSPHERESHELLCCGRADIENT(CHF_FRA(a_gradFab),
                                  CHF_CONST_FRA(diffFab),
                                  CHF_CONST_FRA(xiFab),
                                  CHF_BOX(a_box));
}

//-----------------------------------------------------------------------

void
CubedSphereShellPanelCS::magnitude(
                                   const FArrayBox& a_vecFab,
                                   FArrayBox& a_magFab,
                                   const Box& a_box) const
{
  CH_assert(a_vecFab.box().contains(a_box));
  CH_assert(a_magFab.box().contains(a_box));
  FArrayBox xiFab(a_box, SpaceDim);
  getCenterMappedCoordinates(xiFab, a_box);

  int ncomp = a_magFab.nComp();
  // Order of components in a_vecFab will be:
  // ncomp differences with dimension 0,
  // ncomp differences with dimension 1, etc.
  CH_assert(a_vecFab.nComp() == SpaceDim*ncomp);
  FORT_CUBEDSPHERESHELLMAGNITUDE(CHF_FRA(a_magFab),
                                 CHF_CONST_FRA(a_vecFab),
                                 CHF_CONST_FRA(xiFab),
                                 CHF_BOX(a_box));
}

//-----------------------------------------------------------------------

void
CubedSphereShellPanelCS::magnitudeNatural(
                                          const FArrayBox& a_vecFab,
                                          FArrayBox& a_magFab,
                                          const Box& a_box) const
{
  CH_assert(a_vecFab.nComp() == SpaceDim);
  CH_assert(a_magFab.nComp() == 1);

  CH_assert(a_vecFab.box().contains(a_box));
  CH_assert(a_magFab.box().contains(a_box));
  FArrayBox xiFab(a_box, SpaceDim);
  getCenterMappedCoordinates(xiFab, a_box);

  FArrayBox vecCenUnitFab(a_box, SpaceDim);
  // Convert vector components from natural basis to unit basis.
  FORT_CUBEDSPHERESHELLVECTORTOUNITBASIS(CHF_CONST_FRA(xiFab),
                                         CHF_CONST_FRA(a_vecFab),
                                         CHF_FRA(vecCenUnitFab),
                                         CHF_BOX(a_box));

  FORT_CUBEDSPHERESHELLMAGNITUDE(CHF_FRA(a_magFab),
                                 CHF_CONST_FRA(vecCenUnitFab),
                                 CHF_CONST_FRA(xiFab),
                                 CHF_BOX(a_box));
}

//-----------------------------------------------------------------------

void CubedSphereShellPanelCS::getContravariantMetricJAverage(FArrayBox& a_matrixFab,const int a_dir) const
{
  CH_TIME("CubedSphereShellPanelCS::getContravariantMetricJAverage");

  CH_assert(a_matrixFab.nComp() == SpaceDim*SpaceDim);
  a_matrixFab.setVal(0);
  Box faceBox=a_matrixFab.box();
  FArrayBox xiFab(faceBox,SpaceDim);
  getCenterMappedCoordinates(xiFab,faceBox);
  FORT_CUBEDSPHEREMETRICJPRODAVG(CHF_FRA(a_matrixFab),
                                 CHF_CONST_FRA(xiFab),
                                 CHF_BOX(faceBox),
                                 CHF_CONST_REALVECT(m_dx),
                                 CHF_CONST_REAL(m_height),
                                 CHF_CONST_REAL(m_radius),
                                 CHF_CONST_INT(a_dir));
}

#include "NamespaceFooter.H"
