#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TWISTEDCS_H_
#define _TWISTEDCS_H_

#include "NewFourthOrderCoordSys.H"


#include "NamespaceHeader.H"


/// Calhoun-Leveque "twisted" Coordinate mapping (constant Jacobian)
///
/**
*/
class TwistedCS : public NewFourthOrderCoordSys
{
public:
  /// constructor
  TwistedCS(const RealVect& a_dX,
            const Real& a_R,
            const Real& a_twist);

  /**
     Destructor.
   */
  virtual ~TwistedCS();

  /// given coordinate in mapped space, return its location in real space
  RealVect realCoord(const RealVect& a_Xi) const;

  /// given coordinate in real space, return its location in the mapped space
  RealVect mappedCoord(const RealVect& a_x) const;

  // Calculate the derivative of each coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;


protected:

  // radius for n-sphere of twisted grid
  Real m_R;

  // Amount of twist in the grid
  Real m_theta;

  // convenience
  Real m_Pi;
  Real m_scale;

};


/// factory for TwistedCS
/**
 */
class
TwistedCSFactory : public NewCoordSysFactory
{

public:
  /// constructor
  TwistedCSFactory(const Real& a_R,
                   const Real& a_twist);

  ///
  virtual NewCoordSys* getCoordSys(const ProblemDomain& a_domain,
                                   const RealVect& a_dx) const;

protected:

  ///
  Real m_R;
  Real m_twist;
};

#include "NamespaceFooter.H"

#endif
