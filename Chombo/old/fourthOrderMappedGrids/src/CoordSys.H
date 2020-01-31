#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _COORDSYS_H_
#define _COORDSYS_H_

#include <utility>
using std::pair;

#include "REAL.H"
#include "DisjointBoxLayout.H"
#include "RealVect.H"
#include "IntVect.H"
#include "LevelData.H"
#include "Vector.H"

#include "NamespaceHeader.H"

template <class Tcell,class Tface>
class CoordSys;



/// Virtual base class encapsulating mapped-grid coordinate systems
/**
   The CoordSys class is a virtual base class providing an API
   for performing operations in a mapped grid coordinate space.
   Implementations of specific coordinate systems will be provided by
   derived classes.
*/
template <class T, class Tface>
class CoordSys
{
public:
  /// default constructor
  CoordSys()
  {
  }

  /**
     Destructor.
   */
  virtual ~CoordSys()
  {;}

  ///
  virtual void define(const DisjointBoxLayout& a_grids,
                      const ProblemDomain& a_domain,
                      const RealVect& a_cellSpacing,
                      const IntVect& a_ghostVect) = 0;

  ///
  virtual void regrid(const DisjointBoxLayout& a_newGrids) = 0;

  /// given coordinate in mapped space, return its location in real space
  virtual RealVect realCoord(const RealVect& a_X) const = 0;

  /// given coordinate in real space, return its location in the mapped space
  virtual RealVect mappedCoord(const RealVect& a_x) const = 0;

  /// given coordinates in mapped space, return locations in real space
  virtual void realCoord(T& a_x, const T& a_X) const = 0;

  /// given coordinate in real space, return its location in the mapped space
  virtual void mappedCoord(T& a_X, const T& a_x) const = 0;

  /// Jacobian evaluated at location X in real space
  virtual Real pointwiseJ(const RealVect& a_X) const = 0;

  /// Jacobian evaluated at locations Xi in mapped space
  virtual void pointwiseJ(T& a_J, const T& a_Xi) const = 0;

  /// returns integral of divergence over mapped-grid cells
  virtual void mappedGridDivergence(LevelData<T>& a_divF,
                                    const LevelData<Tface>& a_F) const = 0;

  /// reuturns reference to face-averaged N metric terms
  virtual const LevelData<Tface>& getFaceMetricTerms() const = 0;

  // returns reference to cell volumes
  virtual const LevelData<T>& getCellVolumes() const = 0;

  /// returns second-order averaged cell-centered J
  virtual const LevelData<T>& getJ() const = 0;

  /// returns averaged face-centered 1/J
  virtual const LevelData<Tface>& getJInverse() const =0;

  /// returns averaged face-centered N/J
  virtual const LevelData<Tface>& getNJinverse() const = 0;
};


/// pure-virtual base class defining factory interface for CoordSys
template <class T, class Tface>
class
CoordSysFactory
{


public:

  CoordSysFactory()
  {
  }

  virtual ~CoordSysFactory()
  {
  }

  virtual CoordSys<T,Tface>* getCoordSys(const DisjointBoxLayout& grids,
                                         const ProblemDomain& levelDomain,
                                         const IntVect& a_ghostVect ) const = 0;



};

#include "NamespaceFooter.H"

#endif
