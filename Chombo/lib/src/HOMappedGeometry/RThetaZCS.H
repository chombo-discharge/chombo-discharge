#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _RTHETAZCS_H_
#define _RTHETAZCS_H_

#include "NewFourthOrderCoordSys.H"


#include "NamespaceHeader.H"

/// Basic (possibly stretched) Cartesian fourth-order oordinate system
/**
   The CartesianCS class implements the NewFourthOrderCoordSys class
   for Cartesian mappings
*/
class RThetaZCS : public NewFourthOrderCoordSys
{
public:
  /// constructor
  RThetaZCS(const RealVect& a_dX, const RealVect a_stretch = RealVect::Unit,
            const RealVect& m_origin = RealVect::Zero);

  /**
     Destructor.
   */
  virtual ~RThetaZCS();

  /// given coordinate in mapped space, return its location in real space
  RealVect realCoord(const RealVect& a_Xi) const;

  /// given coordinate in real space, return its location in the mapped space
  RealVect mappedCoord(const RealVect& a_x) const;

  // Calculate the derivative of each coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;

protected:
  // allows for stretched Cartesian grids
  RealVect m_stretch;

  // allows for offsets (allows avoiding singularity at center)
  RealVect m_origin;

  // convenience. It's always good to have some pi around
  Real m_Pi;

};


/// factory for RThetaZCS
/**
 */
class
RThetaZCSFactory : public NewCoordSysFactory
{

public:
  /// constructor
  RThetaZCSFactory(const RealVect& a_stretch = RealVect::Unit,
                   const RealVect& a_origin = RealVect::Zero);

  ///
  virtual NewCoordSys* getCoordSys(const ProblemDomain& a_domain,
                                   const RealVect& a_dx) const;

protected:

  ///
  RealVect m_stretch;
  ///
  RealVect m_origin;
};

#include "NamespaceFooter.H"

#endif
