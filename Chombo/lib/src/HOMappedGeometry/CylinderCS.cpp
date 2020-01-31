#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CylinderCS.H"
#include "CylinderSettingsF_F.H"

#include "NamespaceHeader.H"

CylinderCS::CylinderCS()
{
  m_coordSysVect.resize(NUMBLOCKS, NULL);

  m_mappingBlocks.resize(NUMBLOCKS);
}

CylinderCS::~CylinderCS()
{
}

void
CylinderCS::define(const ProblemDomain& a_levelDomain,
                   const RealVect& a_dx)
{
  const Box& levelBox = a_levelDomain.domainBox();
  int baseLength = levelBox.size(0) / 5;
  CH_assert(levelBox.size(0) == 5 * baseLength);
  CH_assert(levelBox.size(1) == 5 * baseLength);

  IntVect baseLo = IntVect(D_DECL6(0, 0, levelBox.smallEnd(2), 0, 0, 0));
  IntVect baseHi = IntVect(D_DECL6(baseLength-1, baseLength-1, levelBox.bigEnd(2), 0, 0, 0));
  Box baseBox = Box(baseLo, baseHi);

  RealVect baseRealLength = a_dx * Real(baseLength);
  RealVect baseRealLo = RealVect::Zero;
  RealVect baseRealHi = baseRealLength;
  RealVect basis0 = BASISREALV(0);
  RealVect basis1 = BASISREALV(1);

  Vector<IntVect> offsets(NUMBLOCKS);
  offsets[CENTRAL] = IntVect::Zero;
  offsets[XPOS] = 2 * BASISV(0);
  offsets[YPOS] = 2 * BASISV(1);
  offsets[XNEG] = -2 * BASISV(0);
  offsets[YNEG] = -2 * BASISV(1);

  for (int iblock = 0; iblock < NUMBLOCKS; iblock++)
    {
      m_mappingBlocks[iblock] = baseBox + baseLength * offsets[iblock];
    }
  m_gotMappingBlocks = true;

  /*
    Define block boundaries.
   */
  defineBoundaries();

  initializeBlockTransformations();
  // m_coordSysVect[0] = new CylinderBlockCS(...)
}

void
CylinderCS::defineBoundaries()
{
  // Start by setting all block boundaries to be of type BOUNDARY.
  setAllBoundaries(BlockBoundary::BOUNDARY);

  int xLo = 0;
  int xHi = SpaceDim;
  int yLo = 1;
  int yHi = 1 + SpaceDim;
  IntVect clockwise = IntVect(D_DECL6(-1, 1, 1, 1, 1, 1));
  IntVect anticlockwise = IntVect(D_DECL6(1, -1, 1, 1, 1, 1));

  // Block CENTRAL = 0:  [0:1, 0:1]
  setBoundaryFromFaces(CENTRAL, xHi, XPOS, xLo);
  setBoundaryFromFaces(CENTRAL, yHi, YPOS, yLo);
  setBoundaryFromFaces(CENTRAL, xLo, XNEG, xHi);
  setBoundaryFromFaces(CENTRAL, yLo, YNEG, yHi);

  // Block XPOS = 1:  [2:3, 0:1]
  setBoundaryFromFaces(XPOS, xLo, CENTRAL, xHi);
  setBoundaryFromFaces(XPOS, yHi, YPOS, xHi, anticlockwise);
  setBoundaryFromFaces(XPOS, yLo, YNEG, xHi, clockwise);

  // Block YPOS = 2:  [0:1, 2:3]
  setBoundaryFromFaces(YPOS, yLo, CENTRAL, yHi);
  setBoundaryFromFaces(YPOS, xHi, XPOS, yHi, clockwise);
  setBoundaryFromFaces(YPOS, xLo, XNEG, yHi, anticlockwise);

  // Block XNEG = 3:  [-2:-1, 0:1]
  setBoundaryFromFaces(XNEG, xHi, CENTRAL, xLo);
  setBoundaryFromFaces(XNEG, yHi, YPOS, xLo, clockwise);
  setBoundaryFromFaces(XNEG, yLo, YNEG, xLo, anticlockwise);

  // Block YNEG = 4:  [0:1, -2:-1]
  setBoundaryFromFaces(YNEG, yHi, CENTRAL, yLo);
  setBoundaryFromFaces(YNEG, xHi, XPOS, yLo, anticlockwise);
  setBoundaryFromFaces(YNEG, xLo, XNEG, yLo, clockwise);

  m_gotBoundaries = true;
}


void
CylinderCS::blockRemapping(RealVect& a_xi_valid,
                           int& a_n_valid,
                           const RealVect& a_xiSrc,
                           int a_nSrc) const
{
  RealVect X = m_coordSysVect[a_nSrc]->realCoord(a_xiSrc);

  // Now figure out a_n_valid.
  // Must have m_centralRectLo and m_centralRectHi,
  // which are defined in setAllPhysical().
  if (m_centralRectLo <= X && X <= m_centralRectHi) // all components <=
    {
      a_n_valid = 0;
    }
  else if (X[0] >= abs(X[1]))
    {
      a_n_valid = 1;
    }
  else if (X[1] >= abs(X[0]))
    {
      a_n_valid = 2;
    }
  else if (X[0] <= -abs(X[1]))
    {
      a_n_valid = 3;
    }
  else if (X[1] <= -abs(X[0]))
    {
      a_n_valid = 4;
    }

  a_xi_valid = m_coordSysVect[a_n_valid]->mappedCoord(X);
}

void
CylinderCS::setAllPhysical(const RealVect&  a_centerPoint,
                           const RealVect&  a_centralRectSize,
                           const Real&      a_outerRadius)
{
  m_centerPoint = a_centerPoint;
  m_centralRectSize = a_centralRectSize;
  m_outerRadius = a_outerRadius;

  m_centralRectLo = m_centerPoint - m_centralRectSize / 2.;
  m_centralRectHi = m_centerPoint + m_centralRectSize / 2.;

  setFortranCommon();
}

void
CylinderCS::setFortranCommon()
{
  FORT_CYLINDERSETF(CHF_CONST_REALVECT(m_centerPoint),
                    CHF_CONST_REALVECT(m_centralRectSize),
                    CHF_CONST_REAL(m_outerRadius));
}

void
CylinderCSFactory::setAllPhysical(CylinderCS* a_coordSysPtr) const
{
  // CH_assert(m_isDefined);
  a_coordSysPtr->setAllPhysical(m_centerPoint,
                                m_centralRectSize,
                                m_outerRadius);
}

#include "NamespaceFooter.H"
