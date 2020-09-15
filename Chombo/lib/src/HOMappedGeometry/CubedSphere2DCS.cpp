#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CubedSphere2DCS.H"
#include "CubedSphere2DPanelCS.H"
#include "CubedSphereSettingsF_F.H"
#include "CubedSphere2DF_F.H"
#include "MBMiscUtil.H"

#include <cfloat>

#include "NamespaceHeader.H"
// #include "UsingNamespace.H"

///////////////////////////////////////////////////////////////////////////////

CubedSphere2DCS::CubedSphere2DCS()
{
  m_useScriptN = false;
}

CubedSphere2DCS::CubedSphere2DCS(int a_nResolution,
                                 int a_nLayers)
{
  define(a_nResolution, a_nLayers);
  m_useScriptN = false;
}

CubedSphere2DCS::~CubedSphere2DCS()
{
  if (m_gotCoordSysVect)
    {
      for (int i=0; i<6; i++)
        {
          delete m_coordSysVect[i];
        }
    }
}

void
CubedSphere2DCS::define(int a_nResolution,
                        int a_nLayers)
{
  m_realDim = 3; // changed to 2 by setFlatMap(true) if SpaceDim == 2
  m_flatMap = false;

  // set coordSys for each block before heading off to the big define
  m_coordSysVect.resize(6, NULL);
  m_mappingBlocks.resize(6);

  m_nResolution = a_nResolution;
  m_nLayers = (a_nLayers == 0) ? m_nResolution : a_nLayers;

  // Initialize blocks
  RealVect dx = (0.5*M_PI / static_cast<Real>(m_nResolution)) * RealVect::Unit;
  if (SpaceDim == 3)
    {
      dx[2] = 1. / static_cast<Real>(m_nLayers);
    }

  IntVect small = IntVect::Zero;

  for (int i=0; i<6; i++)
  {
    small[0] = i * 2 * m_nResolution;

    CubedSphere2DPanelCS* thisBlockPtr = new CubedSphere2DPanelCS(i, dx, small);
    m_coordSysVect[i] = reinterpret_cast<NewCoordSys*>(thisBlockPtr);

    IntVect big = small + (m_nResolution - 1) * IntVect::Unit;

    m_mappingBlocks[i] = Box(small, big);
  }
  m_gotCoordSysVect = true;
  m_gotMappingBlocks = true;

  // Define boundaries
  defineBoundaries();

  // Initialize block transformations
  initializeBlockTransformations();
}

void
CubedSphere2DCS::define(const ProblemDomain& a_levelDomain,
                        const RealVect& a_dx)
{
  // Set resolution per block
  const Box& levelDomainBox = a_levelDomain.domainBox();
  if (SpaceDim == 2)
    define(levelDomainBox.size(1));
  else // (SpaceDim == 3)
    define(levelDomainBox.size(1), levelDomainBox.size(2));
}

void
CubedSphere2DCS::setFlatMap(bool a_flatMap)
{
  m_flatMap = a_flatMap;
  for (int i=0; i<6; i++)
    {
      ((CubedSphere2DPanelCS*) m_coordSysVect[i])->setFlatMap(m_flatMap);
    }
  if (SpaceDim == 2)
    if (m_flatMap)
      m_realDim = 2;
}

void
CubedSphere2DCS::setRadii(Real a_rMin,
                          Real a_rMax)
{
  m_rMin = a_rMin;
  m_rMax = a_rMax;
  setFortranCommon();
}

void
CubedSphere2DCS::setFortranCommon()
{
  FORT_CUBEDSPHERESETF(CHF_CONST_REAL(m_rMin),
                       CHF_CONST_REAL(m_rMax));
}


void
CubedSphere2DCS::defineBoundaries()
{
  CH_assert(gotMappingBlocks());

  // length in dimension 0 of central block (index 0)
  int baseLength = m_mappingBlocks[0].size(0);
  m_boundaries.resize(6);

  // Start by setting all block boundaries to be of type BOUNDARY.
  for (int iblock = 0; iblock < 6; iblock++)
    {
      Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[iblock];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          blockBoundaries[idir].define(0);
          blockBoundaries[idir + SpaceDim].define(0);
        }
    }

  int xLo = 0;
  int xHi = SpaceDim;
  int yLo = 1;
  int yHi = 1 + SpaceDim;
  IntVect xShift = baseLength*BASISV(0);
  IntVect yShift = baseLength*BASISV(1);
  IntVect permSame = IntVect(D_DECL6(0, 1, 2, 0, 0, 0));
  IntVect permDiff = IntVect(D_DECL6(1, 0, 2, 0, 0, 0));
  IntVect clockwise = IntVect(D_DECL6(-1, 1, 1, 0, 0, 0));
  IntVect anticlockwise = IntVect(D_DECL6(1, -1, 1, 0, 0, 0));
  IntVect flip = IntVect(D_DECL6(-1, -1, 1, 0, 0, 0));
  IndicesTransformation it;

  // Equatorial panel 0
  it.defineFromTranslation(xShift);
  m_boundaries[0][xHi].define(it, 1);
  it.defineFromTranslation(8*xShift - yShift);
  m_boundaries[0][yHi].define(it, 4);
  it.defineFromTranslation(7*xShift);
  m_boundaries[0][xLo].define(it, 3);
  it.defineFromTranslation(10*xShift + yShift);
  m_boundaries[0][yLo].define(it, 5);

  // Equatorial panel 1
  it.defineFromTranslation(xShift);
  m_boundaries[1][xHi].define(it, 2);
  it.defineFromPivot(2*xShift + yShift,
                     9*xShift,
                     permDiff, anticlockwise);
  m_boundaries[1][yHi].define(it, 4);
  it.defineFromTranslation(-xShift);
  m_boundaries[1][xLo].define(it, 0);
  it.defineFromPivot(2*xShift,
                     11*xShift + yShift,
                     permDiff, clockwise);
  m_boundaries[1][yLo].define(it, 5);

  // Equatorial panel 2
  it.defineFromTranslation(xShift);
  m_boundaries[2][xHi].define(it, 3);
  it.defineFromPivot(4*xShift + yShift,
                     9*xShift + yShift,
                     permSame, flip);
  m_boundaries[2][yHi].define(it, 4);
  it.defineFromTranslation(-xShift);
  m_boundaries[2][xLo].define(it, 1);
  it.defineFromPivot(4*xShift,
                     11*xShift,
                     permSame, flip);
  m_boundaries[2][yLo].define(it, 5);

  // Equatorial panel 3
  it.defineFromTranslation(-7*xShift);
  m_boundaries[3][xHi].define(it, 0);
  it.defineFromPivot(7*xShift + yShift,
                     8*xShift,
                     permDiff, clockwise);
  m_boundaries[3][yHi].define(it, 4);
  it.defineFromTranslation(-xShift);
  m_boundaries[3][xLo].define(it, 2);
  it.defineFromPivot(7*xShift,
                     10*xShift + yShift,
                     permDiff, anticlockwise);
  m_boundaries[3][yLo].define(it, 5);

  // North polar panel (4)
  it.defineFromPivot(9*xShift,
                     2*xShift + yShift,
                     permDiff, clockwise);
  m_boundaries[4][xHi].define(it, 1);
  it.defineFromPivot(9*xShift + yShift,
                     4*xShift + yShift,
                     permSame, flip);
  m_boundaries[4][yHi].define(it, 2);
  it.defineFromPivot(8*xShift,
                     7*xShift + yShift,
                     permDiff, anticlockwise);
  m_boundaries[4][xLo].define(it, 3);
  it.defineFromTranslation(-8*xShift + yShift);
  m_boundaries[4][yLo].define(it, 0);

  // South polar panel (5)
  it.defineFromPivot(11*xShift + yShift,
                     2*xShift,
                     permDiff, anticlockwise);
  m_boundaries[5][xHi].define(it, 1);
  it.defineFromTranslation(-10*xShift - yShift);
  m_boundaries[5][yHi].define(it, 0);
  it.defineFromPivot(10*xShift + yShift,
                     7*xShift,
                     permDiff, clockwise);
  m_boundaries[5][xLo].define(it, 3);
  it.defineFromPivot(11*xShift,
                     4*xShift,
                     permSame, flip);
  m_boundaries[5][yLo].define(it, 2);

  m_gotBoundaries = true;
}

RealVect
CubedSphere2DCS::realCoord(const RealVect& a_X) const
{
  // Obtain the panel index from the alpha coordinate
  int nPanel = static_cast<int>((a_X[0] + 0.5 * M_PI) / M_PI);

  // Assert that the panel index is consistent
  CH_assert((nPanel >= 0) && (nPanel <= 5));

  // Submit
  return m_coordSysVect[nPanel]->realCoord(a_X);
}

RealVect
CubedSphere2DCS::mappedCoord(const RealVect& a_x) const
{
  // Determine the panel index
  int nPanel;

  if ((fabs(a_x[0]) > fabs(a_x[1])) && (fabs(a_x[0]) > fabs(a_x[2])))
  {
    if (a_x[0] > 0.0)
      nPanel = 0;
    else
      nPanel = 2;
  }
  else if (fabs(a_x[1]) > fabs(a_x[2]))
  {
    if (a_x[1] > 0.0)
      nPanel = 1;
    else
      nPanel = 3;
  }
  else
  {
    if (a_x[2] > 0.0)
      nPanel = 4;
    else
      nPanel = 5;
  }

  // Return the mapped coord of the appropriate panel
  return m_coordSysVect[nPanel]->mappedCoord(a_x);
}

/*
Real
CubedSphere2DCS::dXdXi(const RealVect& a_X, int a_dirX, int a_dirXi) const
{
  Real value = OldMultiBlockCoordSys::dXdXi(a_X, a_dirX, a_dirXi);

  return value;
}
*/

void CubedSphere2DCS::vecTransMatrix_Pt(
    int a_nSourcePanel,
    int a_nDestPanel,
    Real a_dX,
    Real a_dY,
    Real* a_dT
) const
{
  // Called by vecTransMatrix() and also by vectorBlockTransformation(),
  // which deals with the radial direction in 3D.
  Real dX2 = a_dX * a_dX;
  Real dY2 = a_dY * a_dY;
  Real dNorm = 1.0 / (dX2 + dY2);

  // Equatorial panels
  if ((a_nSourcePanel < 4) && (a_nDestPanel < 4))
  {
    if (a_nSourcePanel == ((a_nDestPanel + 1) % 4))
    {
      a_dT[0] = 1.0;
      a_dT[1] = 0.0;
      a_dT[2] = a_dY * (1.0 + dX2) * dNorm;
      a_dT[3] = - a_dX * (1.0 + dY2) * dNorm;
      return;

    }
    else if (a_nSourcePanel == ((a_nDestPanel + 3) % 4))
    {
      a_dT[0] = 1.0;
      a_dT[1] = 0.0;
      a_dT[2] = - a_dY * (1.0 + dX2) * dNorm;
      a_dT[3] = a_dX * (1.0 + dY2) * dNorm;
      return;

    }
    else
    {
      CH_assert(false);
    }

  // North polar panel is the source
  }
  else if (a_nSourcePanel == 4)
  {
    a_dT[0] = - a_dY * (1.0 + dX2) * dNorm;
    a_dT[1] = a_dX * (1.0 + dY2) * dNorm;

    if (a_nDestPanel == 0)
    {
      a_dT[2] = 0.0;
      a_dT[3] = 1.0;
      return;

    }
    else if (a_nDestPanel == 1)
    {
      a_dT[2] = -1.0;
      a_dT[3] = 0.0;
      return;

    }
    else if (a_nDestPanel == 2)
    {
      a_dT[2] = 0.0;
      a_dT[3] = -1.0;
      return;

    }
    else if (a_nDestPanel == 3)
    {
      a_dT[2] = 1.0;
      a_dT[3] = 0.0;
      return;

    }
    else
    {
      CH_assert(false);
    }

  // North polar panel is the destination
  }
  else if (a_nDestPanel == 4)
  {
    if (a_nSourcePanel == 0)
    {
      a_dT[0] = a_dY * (1.0 + dX2) * dNorm;
      a_dT[1] = - a_dX * (1.0 + dY2) * dNorm;
      a_dT[2] = 0.0;
      a_dT[3] = 1.0;
      return;

    }
    else if (a_nSourcePanel == 1)
    {
      a_dT[0] = 0.0;
      a_dT[1] = -1.0;
      a_dT[2] = a_dY * (1.0 + dX2) * dNorm;
      a_dT[3] = - a_dX * (1.0 + dY2) * dNorm;
      return;

    }
    else if (a_nSourcePanel == 2)
    {
      a_dT[0] = - a_dY * (1.0 + dX2) * dNorm;
      a_dT[1] = a_dX * (1.0 + dY2) * dNorm;
      a_dT[2] = 0.0;
      a_dT[3] = -1.0;
      return;

    }
    else if (a_nSourcePanel == 3)
    {
      a_dT[0] = 0.0;
      a_dT[1] = 1.0;
      a_dT[2] = - a_dY * (1.0 + dX2) * dNorm;
      a_dT[3] = a_dX * (1.0 + dY2) * dNorm;
      return;

    }
    else
    {
      CH_assert(false);
    }

  // South polar panel is the source
  }
  else if (a_nSourcePanel == 5)
  {
    a_dT[0] = a_dY * (1.0 + dX2) * dNorm;
    a_dT[1] = - a_dX * (1.0 + dY2) * dNorm;

    if (a_nDestPanel == 0)
    {
      a_dT[2] = 0.0;
      a_dT[3] = 1.0;
      return;

    }
    else if (a_nDestPanel == 1)
    {
      a_dT[2] = 1.0;
      a_dT[3] = 0.0;
      return;

    }
    else if (a_nDestPanel == 2)
    {
      a_dT[2] = 0.0;
      a_dT[3] = -1.0;
      return;

    }
    else if (a_nDestPanel == 3)
    {
      a_dT[2] = -1.0;
      a_dT[3] = 0.0;
      return;

    }
    else
    {
      CH_assert(false);
    }

  // South polar panel is the destination
  }
  else if (a_nDestPanel == 5)
  {
    if (a_nSourcePanel == 0)
    {
      a_dT[0] = - a_dY * (1.0 + dX2) * dNorm;
      a_dT[1] = a_dX * (1.0 + dY2) * dNorm;
      a_dT[2] = 0.0;
      a_dT[3] = 1.0;
      return;

    }
    else if (a_nSourcePanel == 1)
    {
      a_dT[0] = 0.0;
      a_dT[1] = 1.0;
      a_dT[2] = a_dY * (1.0 + dX2) * dNorm;
      a_dT[3] = - a_dX * (1.0 + dY2) * dNorm; // corrected sign, petermc, 2 Aug 2010
      return;

    }
    else if (a_nSourcePanel == 2)
    {
      a_dT[0] = a_dY * (1.0 + dX2) * dNorm;
      a_dT[1] = - a_dX * (1.0 + dY2) * dNorm;
      a_dT[2] = 0.0;
      a_dT[3] = -1.0;
      return;

    }
    else if (a_nSourcePanel == 3)
    {
      a_dT[0] = 0.0;
      a_dT[1] = -1.0;
      a_dT[2] = - a_dY * (1.0 + dX2) * dNorm;
      a_dT[3] = a_dX * (1.0 + dY2) * dNorm;
      return;

    }
    else
    {
      CH_assert(false);
    }

  }
  else
  {
    CH_assert(false);
  }
}


void CubedSphere2DCS::vecTransMatrix(
    int a_nSourcePanel,
    int a_nDestPanel,
    int a_nAlphaIx,
    int a_nBetaIx,
    Real a_dDeltaA,
    Real* a_dAvgT,
    Real* a_dAlphaGradT,
    Real* a_dBetaGradT
) const
{
  // This is for 2D.
  int c;
  int m;
  int n;

  // Centroid coordinates
  Real dCentA = (static_cast<Real>(a_nAlphaIx) + 0.5) * a_dDeltaA - 0.25 * M_PI;
  Real dCentB = (static_cast<Real>(a_nBetaIx)  + 0.5) * a_dDeltaA - 0.25 * M_PI;

  Real dCentX = tan(dCentA);
  Real dCentY = tan(dCentB);

  // Sampled matrix
  Real dPtT[4];

  for (c = 0; c < 4; c++)
  {
    a_dAvgT[c] = 0.0;
  }

  // Calculate the average of the transformation matrix to O(dA^4) using
  // four point Gaussian quadrature.
  for (m = -1; m < 2; m += 2)
  {
  for (n = -1; n < 2; n += 2)
  {
    vecTransMatrix_Pt(
      a_nSourcePanel,
      a_nDestPanel,
      tan(dCentA + 0.5 * a_dDeltaA * static_cast<Real>(m) / sqrt(3.0)),
      tan(dCentB + 0.5 * a_dDeltaA * static_cast<Real>(n) / sqrt(3.0)),
      &(dPtT[0]));

    for (c = 0; c < 4; c++)
    {
      a_dAvgT[c] += 0.25 * dPtT[c];
    }
  }
  }

  // Calculate the gradients of the transformation matrix using centered
  // difference formulas.
  vecTransMatrix_Pt(
    a_nSourcePanel,
    a_nDestPanel,
    tan(dCentA + 0.5 * a_dDeltaA),
    dCentY,
    &(a_dAlphaGradT[0]));

  vecTransMatrix_Pt(
    a_nSourcePanel,
    a_nDestPanel,
    tan(dCentA - 0.5 * a_dDeltaA),
    dCentY,
    &(dPtT[0]));

  for (c = 0; c < 4; c++)
  {
    a_dAlphaGradT[c] = (a_dAlphaGradT[c] - dPtT[c]) / a_dDeltaA;
  }

  vecTransMatrix_Pt(
    a_nSourcePanel,
    a_nDestPanel,
    dCentX,
    tan(dCentB + 0.5 * a_dDeltaA),
    &(a_dBetaGradT[0]));

  vecTransMatrix_Pt(
    a_nSourcePanel,
    a_nDestPanel,
    dCentX,
    tan(dCentB - 0.5 * a_dDeltaA),
    &(dPtT[0]));

  for (c = 0; c < 4; c++)
  {
    a_dBetaGradT[c] = (a_dBetaGradT[c] - dPtT[c]) / a_dDeltaA;
  }
}

VectorTransformation CubedSphere2DCS::vectorBlockTransformation(
    int a_nDst,
    const RealVect& a_xiSrc,
    int a_nSrc
) const
{
  VectorTransformation trans = VectorTransformation::Identity;
  // In 3D, trans is still the identity transformation in dimension 2.
  if (a_nDst != a_nSrc)
    {
      Real dT[4];

      // Calculate the pointwise transform
      vecTransMatrix_Pt(
                        a_nSrc,
                        a_nDst,
                        tan(a_xiSrc[0] - 0.25 * M_PI),
                        tan(a_xiSrc[1] - 0.25 * M_PI),
                        dT);

      // Set the components of the VectorTransformation
      trans.setComponent(0, 0, dT[0]);
      trans.setComponent(0, 1, dT[1]);
      trans.setComponent(1, 0, dT[2]);
      trans.setComponent(1, 1, dT[3]);
    }
  return trans;
}

void CubedSphere2DCS::cartesianFromMapped(
  Real& xx,
  Real& yy,
  Real& zz,
  const RealVect& a_xiSrc,
  int a_nSrc
) const
{
  CH_assert((a_nSrc >= 0) && (a_nSrc <= 5));
  FORT_CUBEDSPHERE2DMAPPEDTOCARTESIAN(CHF_REAL(xx),
                                      CHF_REAL(yy),
                                      CHF_REAL(zz),
                                      CHF_CONST_REALVECT(a_xiSrc),
                                      CHF_CONST_INT(a_nSrc));
  /*
  Real sx, sy, sz;

  Real dX, dY;

  dX = tan(a_xiSrc[0] - 0.25 * M_PI);
  dY = tan(a_xiSrc[1] - 0.25 * M_PI);

  // Convert to relative Cartesian coordinates
  sz = 1.0 / sqrt(1.0 + dX * dX + dY * dY);
  sx = sz * dX;
  sy = sz * dY;

  // Convert to full Cartesian coordinates
  switch(a_nSrc)
  {
    case 0:
      yy = sx;
      zz = sy;
      xx = sz;
      break;

    case 1:
      xx = -sx;
      zz = sy;
      yy = sz;
      break;

    case 2:
      yy = -sx;
      zz = sy;
      xx = -sz;
      break;

    case 3:
      xx = sx;
      zz = sy;
      yy = -sz;
      break;

    case 4:
      yy = sx;
      xx = -sy;
      zz = sz;
      break;

    case 5:
      yy = sx;
      xx = sy;
      zz = -sz;
      break;
  }
#if CH_SPACEDIM == 3
  // Up to now, assumes physical radius is 1.
  Real rReal = radialMappedToReal(a_xiSrc[2]);
  xx *= rReal;
  yy *= rReal;
  zz *= rReal;
#endif
  */
}

void CubedSphere2DCS::mappedFromCartesian(
  RealVect& a_xiDst,
  int a_nDst,
  Real xx,
  Real yy,
  Real zz
) const
{
  CH_assert((a_nDst >= 0) && (a_nDst <= 5));
  FORT_CUBEDSPHERE2DCARTESIANTOMAPPED(CHF_REALVECT(a_xiDst),
                                      CHF_CONST_REAL(xx),
                                      CHF_CONST_REAL(yy),
                                      CHF_CONST_REAL(zz),
                                      CHF_CONST_INT(a_nDst));
  /*
  a_xiDst[0] = (0.25 + a_nDst*1.) * M_PI;
  a_xiDst[1] = 0.25 * M_PI;
  switch(a_nDst)
  {
    case 0:
      a_xiDst[0] += atan(yy/xx);
      a_xiDst[1] += atan(zz/xx);
      break;

    case 1:
      a_xiDst[0] += atan(-xx/yy);
      a_xiDst[1] += atan(zz/yy);
      break;

    case 2:
      a_xiDst[0] += atan(yy/xx);
      a_xiDst[1] += atan(-zz/xx);
      break;

    case 3:
      a_xiDst[0] += atan(-xx/yy);
      a_xiDst[1] += atan(-zz/yy);
      break;

    case 4:
      a_xiDst[0] += atan(yy/zz);
      a_xiDst[1] += atan(-xx/zz);
      break;

    case 5:
      a_xiDst[0] += atan(-yy/zz);
      a_xiDst[1] += atan(-xx/zz);
      break;
  }
#if CH_SPACEDIM == 3
  // Up to now, only ratios of xx, yy, zz have been taken.
  Real rReal = sqrt(xx*xx + yy*yy + zz*zz);
  a_xiDst[2] = radialRealToMapped(rReal);
#endif
  */
}

void CubedSphere2DCS::blockRemapping(
  RealVect& a_xi_valid,
  int& a_n_valid,
  const RealVect& a_xiSrc,
  int a_nSrc
) const
{

  Real xx, yy, zz;
  Real sx, sy, sz;

  // Convert to Cartesian coordinates
  cartesianFromMapped(xx, yy, zz, a_xiSrc, a_nSrc);

  // Check maximality of the x coordinate
  if ((fabs(xx) > fabs(yy)) && (fabs(xx) > fabs(zz)))
  {
    if (xx > 0)
    {
      a_n_valid = 0;
      sx = yy;
      sy = zz;
      sz = xx;

    }
    else
    {
      a_n_valid = 2;
      sx = -yy;
      sy = zz;
      sz = -xx;
    }

  // Check maximality of y coordinate
  }
  else if (fabs(yy) > fabs(zz))
  {
    if (yy > 0)
    {
      a_n_valid = 1;
      sx = -xx;
      sy = zz;
      sz = yy;

    }
    else
    {
      a_n_valid = 3;
      sx = xx;
      sy = zz;
      sz = -yy;
    }

  // Check maximality of z coordinate
  }
  else
  {
    if (zz > 0)
    {
      a_n_valid = 4;
      sx = yy;
      sy = -xx;
      sz = zz;

    }
    else
    {
      a_n_valid = 5;
      sx = yy;
      sy = xx;
      sz = -zz;
    }
  }

  // (DBG) Check the free z coordinate
  CH_assert(sz > 0.0);

  // Use panel information to calculate (alpha, beta) coords
  a_xi_valid[0] = atan(sx / sz)
    + M_PI * static_cast<double>(a_n_valid) + 0.25 * M_PI;
  a_xi_valid[1] = atan(sy / sz) + 0.25 * M_PI;
#if CH_SPACEDIM == 3
  a_xi_valid[2] = a_xiSrc[2];
#endif
}

RealVect CubedSphere2DCS::blockRemappingGeneral(
  int a_nDst,
  const RealVect & a_xiSrc,
  int a_nSrc
) const
{
  Real xx, yy, zz;

  RealVect xiDst;

  // Convert to Cartesian coordinates
  cartesianFromMapped(xx, yy, zz, a_xiSrc, a_nSrc);

  // Convert from Cartesian coordinates
  mappedFromCartesian(xiDst, a_nDst, xx, yy, zz);

  return xiDst;
}

Vector<RealVect>
CubedSphere2DCS::displacements(const Vector<RealVect>&   a_dstCoords,
                               const Vector<int>&        a_dstBlocks,
                               const RealVect&           a_srcCoords,
                               int                       a_srcBlock) const
{ CH_TIME("CubedSphere2DCS::displacements");
  Vector<Real> dstCoordsAll = VectorRealVectToReal(a_dstCoords);
  /*
  int dstNum = a_dstCoords.size();
  Vector<Real> dstCoordsAll(dstNum * SpaceDim);
  int arri = 0;
  for (int i = 0; i < dstNum; i++)
    {
      const RealVect& dstCoordsThis = a_dstCoords[i];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          dstCoordsAll[arri + idir] = dstCoordsThis[idir];
        }
      arri += SpaceDim;
    }
  */
  Vector<Real> dispAll(a_dstCoords.size() * SpaceDim);
  FORT_CUBEDSPHERE2DDISPLACEMENTS(CHF_VR(dispAll),
                                  CHF_CONST_VR(dstCoordsAll),
                                  CHF_CONST_VI(a_dstBlocks),
                                  CHF_CONST_REALVECT(a_srcCoords),
                                  CHF_CONST_INT(a_srcBlock));
  Vector<RealVect> returnVec = VectorRealToRealVect(dispAll);
  /*
  Vector<RealVect> returnVec(dstNum);
  arri = 0;
  for (int i = 0; i < dstNum; i++)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          returnVec[i][idir] = dispAll[arri + idir];
        }
      arri += SpaceDim;
    }
  */
  /*
  Real xSrc, ySrc, zSrc;
  cartesianFromMapped(xSrc, ySrc, zSrc, a_srcCoords, a_srcBlock);
#if CH_SPACEDIM == 3
  Real rRealSrc = sqrt(xSrc*xSrc + ySrc*ySrc + zSrc*zSrc);
  // Get (xSrc, ySrc, zSrc) on unit sphere.
  xSrc /= rRealSrc;
  ySrc /= rRealSrc;
  zSrc /= rRealSrc;
#endif
  
  // First find the matrix that rotates points on the unit sphere
  // in such a way that (xSrc, ySrc, zSrc) gets rotated to (1, 0, 0).

  // Because of a singularity at (-1, 0, 0), if xSrc < 0 then
  // first map (xSrc, ySrc, zSrc) to (-1, 0, 0), then later to (1, 0, 0).
  bool goThruBack = (xSrc < 0);
  Real xFinal = 1.;
  if (goThruBack) xFinal = -1.;
  Real denom = sqrt((xSrc + xFinal)*(xSrc + xFinal) + ySrc*ySrc + zSrc*zSrc);
  // (xAxis, yAxis, zAxis) is unit vector in direction of axis of rotation
  Real xAxis = (xSrc + xFinal)/denom;
  Real yAxis = ySrc / denom;
  Real zAxis = zSrc / denom;
  int len = a_dstCoords.size();
  Vector<RealVect> returnVec(len);
  for (int i = 0; i < len; i++)
    {
      const RealVect& dstXi = a_dstCoords[i];
      int dstBlock = a_dstBlocks[i];
      Real xDst, yDst, zDst;
      cartesianFromMapped(xDst, yDst, zDst, dstXi, dstBlock);
#if CH_SPACEDIM == 3
      Real rRealDst = sqrt(xDst*xDst + yDst*yDst + zDst*zDst);
      // Get (xDst, yDst, zDst) on unit sphere.
      xDst /= rRealDst;
      yDst /= rRealDst;
      zDst /= rRealDst;
#endif
      Real xAxisDst = xAxis * xDst;
      Real yAxisDst = yAxis * yDst;
      Real zAxisDst = zAxis * zDst;
      // Real xRot = (2. * xAxis) * (xAxisDst + yAxisDst + zAxisDst) - xDst;
      Real yRot = (2. * yAxis) * (xAxisDst + yAxisDst + zAxisDst) - yDst;
      if (goThruBack) yRot = -yRot;
      Real zRot = (2. * zAxis) * (xAxisDst + yAxisDst + zAxisDst) - zDst;
      // Convert (xRot, yRot, zRot) to latitude and longitude.
      Real lat = acos(zRot);
      Real lon = asin(yRot  / sin(lat));
      // xi is longitude and latitude shifted.
      RealVect xi;
      xi[0] = lon;
      xi[1] = 0.5 * M_PI - lat;
#if CH_SPACEDIM == 3
      // Now xi[0] and xi[1] are angles,
      // so multiply by a radius to get arc length.
      Real rRealMean = (rRealDst + rRealSrc) * 0.5;
      xi[0] *= rRealMean;
      xi[1] *= rRealMean;
      // Recover radial displacement, as component 2.
      xi[2] = rRealDst - rRealSrc;
#endif
      returnVec[i] = xi;
    }
  */
  return returnVec;
}

void
CubedSphere2DCS::separateVolFlux(LevelData<FluxBox>& a_flux) const
{
  const DisjointBoxLayout& layout = a_flux.disjointBoxLayout();
  DataIterator dit = layout.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box bx = layout[dit];
      int blockNum = whichBlock(bx);
      FluxBox& flub = a_flux[dit];
      switch (blockNum)
        {
        case 1:
        case 3:
        case 5:
          flub[0].negate();
        }
      if (blockNum > 1)
        flub[1].negate();
    }
}

Real
CubedSphere2DCS::radialMappedToReal(Real a_rMapped) const
{
  Real rReal;
  FORT_SUBRMAPPEDTOREAL(CHF_REAL(rReal),
                        CHF_CONST_REAL(a_rMapped));
  return rReal;
}

Real
CubedSphere2DCS::radialRealToMapped(Real a_rReal) const
{
  Real rMapped;
  FORT_SUBRREALTOMAPPED(CHF_REAL(rMapped),
                        CHF_CONST_REAL(a_rReal));
  return rMapped;
}


// -- begin factory implementations ---------------------------

MultiBlockCoordSys*
CubedSphere2DCSFactory::getCoordSys(const ProblemDomain& a_levelDomain,
                                    const RealVect& a_dx) const
{
  CubedSphere2DCS* coordSysPtr = new CubedSphere2DCS();
  coordSysPtr->define(a_levelDomain, a_dx);
  coordSysPtr->setFlatMap(m_flatMap);
  if (SpaceDim >= 3)
    {
      coordSysPtr->setRadii(m_rMin, m_rMax);
    }
  return (static_cast <MultiBlockCoordSys*> (coordSysPtr));
}

#include "NamespaceFooter.H"
