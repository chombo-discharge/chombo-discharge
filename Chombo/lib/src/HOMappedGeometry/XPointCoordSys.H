#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _XPOINTCOORDSYS_H_
#define _XPOINTCOORDSYS_H_

#define XPOINT_BLOCK_SEPARATION 12

#include "MultiBlockCoordSys.H"
#include "NamespaceHeader.H"

class XPointCoordSys : public MultiBlockCoordSys
{
public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor.
  */
  XPointCoordSys();

  XPointCoordSys(const ProblemDomain& a_levelDomain,
                 const RealVect& a_dx);

  void define(const ProblemDomain& a_levelDomain,
              const RealVect& a_dx);

  ///
  /** Destructor.
   */
  virtual ~XPointCoordSys();

  /*@}*/

  // virtual void define( Vector<MultiBlockCoordSys *>& geom_vec);

  /// block mapping conversion function
  /** given xi_k in block k, returns which block (k_valid)
      contains that location in its valid region and its       
      location in that block's computational coordinate space.
  */
   virtual void blockRemapping(RealVect& xi_valid,
                               int& n_valid,
                               const RealVect& xiSrc,
                               int nSrc) const;

   virtual void defineBoundaries();

  // virtual IntVect getDecomposition(int block) const {return m_decomp;}

protected:

  // IntVect m_decomp;

   string m_geom_type;
};


class XPointCoordSysFactory : public MultiBlockCoordSysFactory
{
public:

  XPointCoordSysFactory()
  {
  }

  virtual ~XPointCoordSysFactory()
  {
  }

  virtual MultiBlockCoordSys* getCoordSys(const ProblemDomain& levelDomain,
                                          const RealVect& a_dx) const;

protected:

};

#include "NamespaceFooter.H"

#endif
