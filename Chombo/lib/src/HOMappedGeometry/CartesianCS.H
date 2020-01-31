#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CARTESIANCS_H_
#define _CARTESIANCS_H_

#include "NewFourthOrderCoordSys.H"


#include "NamespaceHeader.H"

/// Basic (possibly stretched) Cartesian fourth-order oordinate system
/**
   The CartesianCS class implements the NewFourthOrderCoordSys class
   for Cartesian mappings
*/
class CartesianCS : public NewFourthOrderCoordSys
{
public:
  CartesianCS()
  {
  }

  /// constructor
  CartesianCS(const RealVect& a_dX,
              const RealVect& a_origin = RealVect::Zero,
              const RealVect& a_stretch = RealVect::Unit);

  /**
     Destructor.
   */
  virtual ~CartesianCS();

  /// given coordinate in mapped space, return its location in real space
  RealVect realCoord(const RealVect& a_X) const;

  /// given coordinate in real space, return its location in the mapped space
  RealVect mappedCoord(const RealVect& a_x) const;

  // Calculate the derivative of each coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;

protected:
  // Allow for translations
  RealVect m_origin;
  // allows for stretched Cartesian grids
  RealVect m_stretch;

};


/// factory for CartesianCS
/**
 */
class
CartesianCSFactory : public NewCoordSysFactory
{

public:
  /// constructor
  CartesianCSFactory(const RealVect& a_origin = RealVect::Zero,
                     const RealVect& a_stretch = RealVect::Unit);

  ///
  virtual NewCoordSys* getCoordSys(const ProblemDomain& a_domain,
                                   const RealVect& a_dx) const;

protected:

  ///
  RealVect m_origin;
  ///
  RealVect m_stretch;
};


#include "NamespaceFooter.H"

#endif
