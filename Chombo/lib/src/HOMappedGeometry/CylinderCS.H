#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// CylinderCS.H
// petermc, 22 Apr 2010

#ifndef _CYLINDERCS_H_
#define _CYLINDERCS_H_

#include "MultiBlockCoordSys.H"
#include "NamespaceHeader.H"

/// MultiBlockCoordSys with a central cubic box and four sectors.
class CylinderCS : public MultiBlockCoordSys
/** MultiBlockCoordSys with a central cubic box and four sectors.
 */
{
public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor.
  */
  CylinderCS();

  ///
  /** Destructor.
   */
  virtual ~CylinderCS();

  virtual void define(const ProblemDomain& a_levelDomain,
                      const RealVect& a_dx);

  /*@}*/

  /// block mapping conversion function
  /** given xi_k in block k, returns which block (k_valid)
      contains that location in its valid region and its
      location in that block's computational coordinate space.
  */
  virtual void blockRemapping(RealVect& a_xi_valid,
                              int& a_n_valid,
                              const RealVect& a_xiSrc,
                              int a_nSrc) const;

  void setAllPhysical(const RealVect&  a_centerPoint,
                      const RealVect&  a_centralRectSize,
                      const Real&      a_outerRadius);

  void setCenterPoint(const RealVect& a_centerPoint)
  { m_centerPoint = a_centerPoint; }

  void setCentralRectSize(const RealVect& a_centralRectSize)
  { m_centralRectSize = a_centralRectSize; }

  void setOuterRadius(Real a_outerRadius)
  { m_outerRadius = a_outerRadius; }

protected:

  const static int NUMBLOCKS = 5;

  const static int CENTRAL = 0;
  const static int XPOS = 1;
  const static int YPOS = 2;
  const static int XNEG = 3;
  const static int YNEG = 4;

  RealVect m_centerPoint;

  RealVect m_centralRectSize;

  Real m_outerRadius;

  RealVect m_centralRectLo;

  RealVect m_centralRectHi;

  // call only after m_centerPoint, m_centralRectSize, m_outerRadius
  // have already been set.
  void setFortranCommon();

  // define m_boundaries
  void defineBoundaries();

private:

  void operator = (const CylinderCS&)
  {
  }

  CylinderCS(const CylinderCS&)
  {
  }
};


class CylinderCSFactory : public MultiBlockCoordSysFactory
{

public:

  CylinderCSFactory()
  {
  }

  virtual ~CylinderCSFactory()
  {
  }

  virtual MultiBlockCoordSys* getCoordSys(const ProblemDomain& a_levelDomain,
                                          const RealVect& a_dx) const = 0;

  void setCenterPoint(const RealVect& a_centerPoint)
  { m_centerPoint = a_centerPoint; }

  void setCentralRectSize(const RealVect& a_centralRectSize)
  { m_centralRectSize = a_centralRectSize; }

  void setOuterRadius(Real a_outerRadius)
  { m_outerRadius = a_outerRadius; }

protected:

  RealVect m_centerPoint;

  RealVect m_centralRectSize;

  Real m_outerRadius;

  virtual void setAllPhysical(CylinderCS* a_coordSysPtr) const;
};

#include "NamespaceFooter.H"

#endif
