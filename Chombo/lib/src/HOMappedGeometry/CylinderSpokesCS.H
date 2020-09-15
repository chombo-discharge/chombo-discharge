#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// CylinderSpokesCS.H
// petermc, 8 Oct 2013

#ifndef _CYLINDERSPOKESCS_H_
#define _CYLINDERSPOKESCS_H_

#include "CylinderCS.H"
#include "NamespaceHeader.H"

/// MultiBlockCoordSys with a central cubic box and four sectors.
class CylinderSpokesCS : public CylinderCS
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
  CylinderSpokesCS();

  ///
  /** Destructor.
   */
  virtual ~CylinderSpokesCS();


  virtual void define(const ProblemDomain& a_levelDomain,
                      const RealVect& a_dx);

  /*@}*/

protected:

private:

  void operator=(const CylinderSpokesCS&);

  CylinderSpokesCS(const CylinderSpokesCS&);
};


class CylinderSpokesCSFactory : public CylinderCSFactory
{

public:

  CylinderSpokesCSFactory()
  {
  }

  virtual ~CylinderSpokesCSFactory()
  {
  }

  virtual MultiBlockCoordSys* getCoordSys(const ProblemDomain& levelDomain,
                                          const RealVect& a_dx) const;
};

#include "NamespaceFooter.H"

#endif
