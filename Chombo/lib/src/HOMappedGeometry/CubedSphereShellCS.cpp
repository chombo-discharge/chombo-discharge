#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CubedSphereShellCS.H"
#include "CubedSphereShellPanelCS.H"
#include "CubedSphereShellF_F.H"
#include "MBMiscUtil.H"
#include "CONSTANTS.H"

#include <cfloat>

#include "NamespaceHeader.H"
// #include "UsingNamespace.H"

///////////////////////////////////////////////////////////////////////////////

CubedSphereShellCS::CubedSphereShellCS()
{
  m_useScriptN = false;
}

CubedSphereShellCS::CubedSphereShellCS(int a_nResolution, int a_nLayers, 
    Real a_height, Real a_radius)
{
  define(a_nResolution, a_nLayers, a_height, a_radius);
  m_useScriptN = false;
}

CubedSphereShellCS::~CubedSphereShellCS()
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
CubedSphereShellCS::define(int a_nResolution, int a_nLayers, 
    Real a_height, Real a_radius)
{
  m_realDim = 3; // changed to SpaceDim by setFlatMap(true)
  m_flatMap = false;

  // set coordSys for each block before heading off to the big define
  m_coordSysVect.resize(6, NULL);
  m_mappingBlocks.resize(6);

  m_nResolution = a_nResolution;

  // Initialize blocks
  m_dxVect = (0.5*M_PI / static_cast<Real>(m_nResolution)) * RealVect::Unit;

  CH_assert(a_height > 0);
  m_height = a_height;
  CH_assert(a_radius > 0);
  m_radius = a_radius;

  // Define the vertical grid spacing if we're in 3D
  if (SpaceDim == 3)
  {
    CH_assert(a_nLayers > 0);
    m_nLayers = a_nLayers;
    m_dxVect[SpaceDim-1] = 1.0 / (Real) m_nLayers;
  }

  IntVect small = IntVect::Zero;

  for (int i=0; i<6; i++)
  {
    small[0] = i * 2 * m_nResolution;

    CubedSphereShellPanelCS* thisBlockPtr = 
      new CubedSphereShellPanelCS(i, m_dxVect, small);
    if (SpaceDim == 3)
    {
      thisBlockPtr->setHeight(m_height);
      thisBlockPtr->setRadius(m_radius);
    }
    m_coordSysVect[i] = reinterpret_cast<NewCoordSys*>(thisBlockPtr);

    IntVect big(D_DECL6(m_nResolution-1, m_nResolution-1, m_nLayers-1, 1, 1, 1));
    big += small;

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
CubedSphereShellCS::define(const ProblemDomain& a_levelDomain,
                        const RealVect& a_dx)
{
  // Ignore a_dx, this mapping defines its own
  // Set resolution per block
  const Box& levelDomainBox = a_levelDomain.domainBox();
  if (SpaceDim == 2)
    define(levelDomainBox.size(1));
  else // (SpaceDim == 3)
    define(levelDomainBox.size(1), levelDomainBox.size(2));
}

void
CubedSphereShellCS::setFlatMap(bool a_flatMap)
{
  m_flatMap = a_flatMap;
  if (m_flatMap)
    m_realDim = SpaceDim;
  else
    m_realDim = 3;

  for (int i=0; i<6; i++)
    {
      ((CubedSphereShellPanelCS*) m_coordSysVect[i])->setFlatMap(m_flatMap);
    }
}

void
CubedSphereShellCS::setHeight(Real a_height)
{
  CH_assert(a_height > 0);
  m_height = a_height;
  for (int i=0; i<6; i++)
    {
      ((CubedSphereShellPanelCS*) m_coordSysVect[i])->setHeight(m_height);
    }
}
 
  
Real
CubedSphereShellCS::getHeight() const
{
  CH_assert(m_height > 0);
  return m_height;
}
  

Real
CubedSphereShellCS::getRadius() const
{
  CH_assert(m_radius > 0);
  return m_radius;
}


void
CubedSphereShellCS::setRadius(Real a_radius)
{
  CH_assert(a_radius > 0);
  m_radius = a_radius;
  for (int i=0; i<6; i++)
    {
      ((CubedSphereShellPanelCS*) m_coordSysVect[i])->setRadius(m_radius);
    }
}

void 
CubedSphereShellCS::setVerticalMap(RefCountedPtr<Spline1DMapping> a_map)
{
  m_verticalMap = a_map;
  for (int i=0; i<6; i++)
    ((CubedSphereShellPanelCS*) m_coordSysVect[i])->setVerticalMap(m_verticalMap);
}

void
CubedSphereShellCS::defineBoundaries()
{
  // Start by setting all block boundaries to be of type BOUNDARY.
  setAllBoundaries(BlockBoundary::BOUNDARY);

  // length in dimension 0 of central block (index 0)
  int baseLength = m_mappingBlocks[0].size(0);

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


void 
CubedSphereShellCS::cssXiFromLonLat2D(
    RealVect a_llr, RealVect& a_xi, int& a_nPanel) const
{
  // Default panel to unattainable value
  a_nPanel = -1;

  // Translate from Lon Lat coordinates to XYZ on a unit sphere
  double xx, yy, zz, pm;
  double lonS = a_llr[0]; // - PI/4.0;
  double lat = a_llr[1];
  xx = cos(lonS) * cos(lat);
  yy = sin(lonS) * cos(lat);
  zz = sin(lat);
  pm = Max(fabs(xx), Max(fabs(yy), fabs(zz)));

  // Check the x coordinate
  double sx, sy, sz;
  if (pm == fabs(xx)) {
    if (xx > 0) 
    {
      a_nPanel = 0;
      sx = yy;
      sy = zz;
      sz = xx;
    } 
    else 
    {
      a_nPanel = 2;
      sx = -yy;
      sy = zz;
      sz = -xx;
    }
  }

  // Check maximality of the y coordinate
  if (pm == fabs(yy)) 
  {
    if (yy > 0) 
    {
      a_nPanel = 1;
      sx = -xx;
      sy = zz;
      sz = yy;
    } 
    else 
    {
      a_nPanel = 3;
      sx = xx;
      sy = zz;
      sz = -yy;
    }
  }

  // Check maximality of the z coordinate
  if (pm == fabs(zz)) 
  {
    if (zz > 0) 
    {
      a_nPanel = 4;
      sx = yy;
      sy = -xx;
      sz = zz;
    } 
    else 
    {
      a_nPanel = 5;
      sx = yy;
      sy = xx;
      sz = -zz;
    }
  }

  // Convert to alpha, beta coordinates, with a multiblock panel shift
  a_xi[0] = atan(sx / sz) + PI*(0.25 + (Real) a_nPanel);
  a_xi[1] = atan(sy / sz) + PI/4.0;
  // a_xi[0] = sx / sz + (PI/2.0)*(0.25 + (Real) a_nPanel);
  // a_xi[1] = sy / sz + (PI/4.0);
}


RealVect
CubedSphereShellCS::realCoord(const RealVect& a_X) const
{
  // Obtain the panel index from the alpha coordinate
  int nPanel = static_cast<int>((a_X[0] + 0.5 * M_PI) / M_PI);

  // Assert that the panel index is consistent
  CH_assert((nPanel >= 0) && (nPanel <= 5));

  // Submit
  return m_coordSysVect[nPanel]->realCoord(a_X);
}

RealVect
CubedSphereShellCS::mappedCoord(const RealVect& a_x) const
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

void CubedSphereShellCS::vecTransMatrix_Pt(
    int a_nSourcePanel,
    int a_nDestPanel,
    Real a_dX,
    Real a_dY,
    Real* a_dT
) const
{
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


void CubedSphereShellCS::vecTransMatrix(
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

VectorTransformation CubedSphereShellCS::vectorBlockTransformation(
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

void CubedSphereShellCS::cartesianFromMapped(
  Real& xx,
  Real& yy,
  Real& zz,
  const RealVect& a_xiSrc,
  int a_nSrc
) const
{
  CH_assert((a_nSrc >= 0) && (a_nSrc <= 5));
  FORT_CUBEDSPHERESHELLMAPPEDTOCARTESIAN(CHF_REAL(xx),
                                         CHF_REAL(yy),
                                         CHF_REAL(zz),
                                         CHF_CONST_REALVECT(a_xiSrc),
                                         CHF_CONST_INT(a_nSrc),
                                         CHF_CONST_REAL(m_radius),
                                         CHF_CONST_REAL(m_height));
}

void CubedSphereShellCS::mappedFromCartesian(
  RealVect& a_xiDst,
  int a_nDst,
  Real xx,
  Real yy,
  Real zz
) const
{
  CH_assert((a_nDst >= 0) && (a_nDst <= 5));

  switch(a_nDst)
  {
    case 0:
      a_xiDst[0] = atan(yy/xx) + 0.25 * M_PI;
      a_xiDst[1] = atan(zz/xx) + 0.25 * M_PI;
      break;

    case 1:
      a_xiDst[0] = atan(-xx/yy) + 1.25 * M_PI;
      a_xiDst[1] = atan(zz/yy) + 0.25 * M_PI;
      break;

    case 2:
      a_xiDst[0] = atan(yy/xx) + 2.25 * M_PI;
      a_xiDst[1] = atan(-zz/xx) + 0.25 * M_PI;
      break;

    case 3:
      a_xiDst[0] = atan(-xx/yy) + 3.25 * M_PI;
      a_xiDst[1] = atan(-zz/yy) + 0.25 * M_PI;
      break;

    case 4:
      a_xiDst[0] = atan(yy/zz) + 4.25 * M_PI;
      a_xiDst[1] = atan(-xx/zz) + 0.25 * M_PI;
      break;

    case 5:
      a_xiDst[0] = atan(-yy/zz) + 5.25 * M_PI;
      a_xiDst[1] = atan(-xx/zz) + 0.25 * M_PI;
      break;
  }
}

void CubedSphereShellCS::blockRemapping(
  RealVect& a_xi_valid,
  int& a_n_valid,
  const RealVect& a_xiSrc,
  int a_nSrc
) const
{

  Real xx, yy, zz;
  Real sx, sy, sz;

  // Convert to Cartesian coordinates on sphere with inner radius 1, outer 2
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

RealVect CubedSphereShellCS::blockRemappingGeneral(
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
CubedSphereShellCS::displacements(const Vector<RealVect>&   a_dstCoords,
                               const Vector<int>&        a_dstBlocks,
                               const RealVect&           a_srcCoords,
                               int                       a_srcBlock) const
{
  Vector<Real> dstCoordsAll = VectorRealVectToReal(a_dstCoords);
  Vector<Real> dispAll(a_dstCoords.size() * SpaceDim);
  FORT_CUBEDSPHERESHELLDISPLACEMENTS(CHF_VR(dispAll),
                                     CHF_CONST_VR(dstCoordsAll),
                                     CHF_CONST_VI(a_dstBlocks),
                                     CHF_CONST_REALVECT(a_srcCoords),
                                     CHF_CONST_INT(a_srcBlock));
  Vector<RealVect> returnVec = VectorRealToRealVect(dispAll);
  /*
  Real xSrc, ySrc, zSrc;
  cartesianFromMapped(xSrc, ySrc, zSrc, a_srcCoords, a_srcBlock);
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
      returnVec[i] = xi;
    }
  */
  return returnVec;
}

void
CubedSphereShellCS::separateVolFlux(LevelData<FluxBox>& a_flux) const
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

// -- begin factory implementations ---------------------------

MultiBlockCoordSys*
CubedSphereShellCSFactory::getCoordSys(const ProblemDomain& a_levelDomain,
                                    const RealVect& a_dx) const
{
  CubedSphereShellCS* coordSysPtr = new CubedSphereShellCS();
  coordSysPtr->define(a_levelDomain, a_dx);
  coordSysPtr->setFlatMap(m_flatMap);
  coordSysPtr->setHeight(m_height);
  coordSysPtr->setRadius(m_radius);
  coordSysPtr->setVerticalMap(m_verticalMap);
  // setAllPhysical(coordSysPtr);
  return (static_cast <MultiBlockCoordSys*> (coordSysPtr));
}

#include "NamespaceFooter.H"
