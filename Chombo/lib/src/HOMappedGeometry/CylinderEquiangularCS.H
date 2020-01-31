#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// CylinderEquiangularCS.H
// petermc, 22 Apr 2010

#ifndef _CYLINDEREQUIANGULARCS_H_
#define _CYLINDEREQUIANGULARCS_H_

#include "CylinderCS.H"
#include "NamespaceHeader.H"

/// MultiBlockCoordSys with a central cubic box and four sectors.
class CylinderEquiangularCS : public CylinderCS
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
  CylinderEquiangularCS();

  ///
  /** Destructor.
   */
  virtual ~CylinderEquiangularCS();


  virtual void define(const ProblemDomain& a_levelDomain,
                      const RealVect& a_dx);

  /*@}*/

protected:

private:

  void operator = (const CylinderEquiangularCS&)
  {
  }

  CylinderEquiangularCS(const CylinderEquiangularCS&)
  {
  }
};


class CylinderEquiangularCSFactory : public CylinderCSFactory
{

public:

  CylinderEquiangularCSFactory()
  {
  }

  virtual ~CylinderEquiangularCSFactory()
  {
  }

  virtual MultiBlockCoordSys* getCoordSys(const ProblemDomain& levelDomain,
                                          const RealVect& a_dx) const;
};

#include "NamespaceFooter.H"

#endif
