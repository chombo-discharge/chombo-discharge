#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CubedSphereSolidCS.H"
#include "CubedSphereSolidBlockCS.H"

// #include <cfloat>

#include "NamespaceHeader.H"

///////////////////////////////////////////////////////////////////////////////

CubedSphereSolidCS::CubedSphereSolidCS()
{
}

CubedSphereSolidCS::CubedSphereSolidCS(int a_nResolution, Real a_r0, Real a_r1)
{
  define(a_nResolution, a_r0, a_r1);
}

CubedSphereSolidCS::~CubedSphereSolidCS()
{
  if (gotCoordSysVect())
    {
      for (int iblock = 0; iblock < NUMBLOCKS; iblock++)
        {
          delete m_coordSysVect[iblock];
        }
    }
}

void
CubedSphereSolidCS::define(int a_nResolution, Real a_r0, Real a_r1)
{
  // set coordSys for each block before heading off to the big define
  m_coordSysVect.resize(NUMBLOCKS, NULL);
  m_mappingBlocks.resize(NUMBLOCKS);

  m_nResolution = a_nResolution;

  // Initialize blocks
  m_dxVect = (1.0 / static_cast<Real>(m_nResolution)) * RealVect::Unit;

  CH_assert(a_r0 > 0);
  m_r0 = a_r0;
  CH_assert(a_r1 > 0);
  m_r1 = a_r1;

  IntVect cornerLo = IntVect::Zero;

  Vector<IntVect> offsets(NUMBLOCKS);
  offsets[XPOS] = 2 * BASISV(0);
  offsets[YPOS] = 2 * BASISV(1);
  offsets[XNEG] = -2 * BASISV(0);
  offsets[YNEG] = -2 * BASISV(1);
  offsets[ZPOS] = 2 * BASISV(2);
  offsets[ZNEG] = -2 * BASISV(2);
  offsets[CENTRAL] = IntVect::Zero;

  Box baseBox(IntVect::Zero, (m_nResolution-1) * IntVect::Unit);
  for (int iblock = 0; iblock < NUMBLOCKS; iblock++)
    {
      m_mappingBlocks[iblock] = baseBox + m_nResolution * offsets[iblock];
      m_coordSysVect[iblock] =
        new CubedSphereSolidBlockCS(iblock, m_dxVect, m_mappingBlocks[iblock]);
    }
  m_gotMappingBlocks = true;
  m_gotCoordSysVect = true;

  setR0(m_r0);
  setR1(m_r1);

  // Define boundaries
  defineBoundaries();

  // Initialize block transformations
  initializeBlockTransformations();
}

void
CubedSphereSolidCS::define(const ProblemDomain& a_levelDomain,
                           const RealVect& a_dx)
{
  // Ignore a_dx, this mapping defines its own
  // Set resolution per block
  const Box& levelDomainBox = a_levelDomain.domainBox();
  int baseLength = levelDomainBox.size(0) / 5;
  CH_assert(levelDomainBox.size(0) == 5 * baseLength);
  CH_assert(levelDomainBox.size(1) == 5 * baseLength);
  CH_assert(levelDomainBox.size(2) == 5 * baseLength);
  define(baseLength);
}

void
CubedSphereSolidCS::setR0(Real a_r0)
{
  CH_assert(a_r0 > 0);
  m_r0 = a_r0;
  CH_assert(gotCoordSysVect());
  for (int iblock = 0; iblock < NUMBLOCKS; iblock++)
    {
      ((CubedSphereSolidBlockCS*) m_coordSysVect[iblock])->setR0(m_r0);
    }
}
 
  
void
CubedSphereSolidCS::setR1(Real a_r1)
{
  CH_assert(a_r1 > 0);
  m_r1 = a_r1;
  CH_assert(gotCoordSysVect());
  for (int iblock = 0; iblock < NUMBLOCKS; iblock++)
    {
      ((CubedSphereSolidBlockCS*) m_coordSysVect[iblock])->setR1(m_r1);
    }
}
 
  
void
CubedSphereSolidCS::defineBoundaries()
{
  // Initialize all boundaries to be physical ones.
  setAllBoundaries(BlockBoundary::BOUNDARY);

  int xLo = 0;
  int xHi = SpaceDim;
  int yLo = 1;
  int yHi = 1 + SpaceDim;
  int zLo = 2;
  int zHi = 2 + SpaceDim;
  IntVect clockwiseXY = IntVect(D_DECL6(-1, 1, 1, 1, 1, 1));
  IntVect anticlockwiseXY = IntVect(D_DECL6(1, -1, 1, 1, 1, 1));
  IntVect clockwiseXZ = IntVect(D_DECL6(-1, 1, 1, 1, 1, 1));
  IntVect anticlockwiseXZ = IntVect(D_DECL6(1, 1, -1, 1, 1, 1));
  IntVect clockwiseYZ = IntVect(D_DECL6(1, -1, 1, 1, 1, 1));
  IntVect anticlockwiseYZ = IntVect(D_DECL6(1, 1, -1, 1, 1, 1));

  // Equatorial block XPOS = 0:  [2:3, 0:1, 0:1]
  setBoundaryFromFaces(XPOS, xLo, CENTRAL, xHi);
  setBoundaryFromFaces(XPOS, yHi, YPOS, xHi, anticlockwiseXY);
  setBoundaryFromFaces(XPOS, yLo, YNEG, xHi, clockwiseXY);
  setBoundaryFromFaces(XPOS, zHi, ZPOS, xHi, anticlockwiseXZ);
  setBoundaryFromFaces(XPOS, zLo, ZNEG, xHi, clockwiseXZ);

  // Equatorial block YPOS = 1:  [0:1, 2:3, 0:1]
  setBoundaryFromFaces(YPOS, yLo, CENTRAL, yHi);
  setBoundaryFromFaces(YPOS, xHi, XPOS, yHi, clockwiseXY);
  setBoundaryFromFaces(YPOS, xLo, XNEG, yHi, anticlockwiseXY);
  setBoundaryFromFaces(YPOS, zHi, ZPOS, yHi, anticlockwiseYZ);
  setBoundaryFromFaces(YPOS, zLo, ZNEG, yHi, clockwiseYZ);

  // Equatorial block XNEG = 2:  [-2:-1, 0:1, 0:1]
  setBoundaryFromFaces(XNEG, xHi, CENTRAL, xLo);
  setBoundaryFromFaces(XNEG, yHi, YPOS, xLo, clockwiseXY);
  setBoundaryFromFaces(XNEG, yLo, YNEG, xLo, anticlockwiseXY);
  setBoundaryFromFaces(XNEG, zHi, ZPOS, xLo, clockwiseXZ);
  setBoundaryFromFaces(XNEG, zLo, ZNEG, xLo, anticlockwiseXZ);

  // Equatorial block YNEG = 3:  [0:1, -2:-1, 0:1]
  setBoundaryFromFaces(YNEG, yHi, CENTRAL, yLo);
  setBoundaryFromFaces(YNEG, xHi, XPOS, yLo, anticlockwiseXY);
  setBoundaryFromFaces(YNEG, xLo, XNEG, yLo, clockwiseXY);
  setBoundaryFromFaces(YNEG, zHi, ZPOS, yLo, clockwiseYZ);
  setBoundaryFromFaces(YNEG, zLo, ZNEG, yLo, anticlockwiseYZ);

  // North polar block ZPOS = 4:  [0:1, 0:1, 2:3]
  setBoundaryFromFaces(ZPOS, zLo, CENTRAL, zHi);
  setBoundaryFromFaces(ZPOS, xHi, XPOS, zHi, clockwiseXZ);
  setBoundaryFromFaces(ZPOS, xLo, XNEG, zHi, anticlockwiseXZ);
  setBoundaryFromFaces(ZPOS, yHi, YPOS, zHi, clockwiseYZ);
  setBoundaryFromFaces(ZPOS, yLo, YNEG, zHi, anticlockwiseYZ);

  // South polar block ZNEG = 5:  [0:1, 0:1, -2:-1]
  setBoundaryFromFaces(ZNEG, zHi, CENTRAL, zLo);
  setBoundaryFromFaces(ZNEG, xHi, XPOS, zLo, anticlockwiseXZ);
  setBoundaryFromFaces(ZNEG, xLo, XNEG, zLo, clockwiseXZ);
  setBoundaryFromFaces(ZNEG, yHi, YPOS, zLo, anticlockwiseYZ);
  setBoundaryFromFaces(ZNEG, yLo, YNEG, zLo, clockwiseYZ);

  // Central cubic block CENTRAL = 6:  [0:1, 0:1, 0:1]
  setBoundaryFromFaces(CENTRAL, xHi, XPOS, xLo);
  setBoundaryFromFaces(CENTRAL, yHi, YPOS, yLo);
  setBoundaryFromFaces(CENTRAL, xLo, XNEG, xHi);
  setBoundaryFromFaces(CENTRAL, yLo, YNEG, yHi);
  setBoundaryFromFaces(CENTRAL, zHi, ZPOS, zLo);
  setBoundaryFromFaces(CENTRAL, zLo, ZNEG, zHi);

  m_gotBoundaries = true;
}

void CubedSphereSolidCS::blockRemapping(
                                        RealVect& a_xi_valid,
                                        int& a_n_valid,
                                        const RealVect& a_xiSrc,
                                        int a_nSrc
                                        ) const
{
  CH_assert(gotCoordSysVect());
  RealVect realCoords = m_coordSysVect[a_nSrc]->realCoord(a_xiSrc);

  // Set indMax to index of component with maximum abs(realCoords).
  int indMax = realCoords.maxDir(true);
  Real valMax = realCoords[indMax];
  if (abs(valMax) <= m_r0)
    {
      // All components of realCoords are in range [-m_r0 : m_r0].
      a_n_valid = CENTRAL;
    }
  else
    {
      // If we get this far, then abs(valMax) > m_r0.
      switch(indMax)
        {
        case 0 :
          {
            if (valMax > 0.)
              a_n_valid = XPOS;
            else
              a_n_valid = XNEG;
            break;
          }
        case 1 :
          {
            if (valMax > 0.)
              a_n_valid = YPOS;
            else
              a_n_valid = YNEG;
            break;
          }
        case 2 :
          {
            if (valMax > 0.)
              a_n_valid = ZPOS;
            else
              a_n_valid = ZNEG;
            break;
          }
        default :
          {
            MayDay::Error("Error in CubedSphereSolidCS::blockRemapping");
          }
        }
    }

  a_xi_valid = m_coordSysVect[a_n_valid]->mappedCoord(realCoords);
}

// -- begin factory implementations ---------------------------

MultiBlockCoordSys*
CubedSphereSolidCSFactory::getCoordSys(const ProblemDomain& a_levelDomain,
                                       const RealVect& a_dx) const
{
  CubedSphereSolidCS* coordSysPtr = new CubedSphereSolidCS();
  coordSysPtr->define(a_levelDomain, a_dx);
  coordSysPtr->setR0(m_r0);
  coordSysPtr->setR1(m_r1);
  // setAllPhysical(coordSysPtr);
  return (static_cast <MultiBlockCoordSys*> (coordSysPtr));
}

#include "NamespaceFooter.H"
