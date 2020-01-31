#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// InterpOnFace.H
// petermc, 20 Aug 2004

#ifndef _INTERPONFACE_H_
#define _INTERPONFACE_H_

#include <cstdlib>
#include <iostream>
#include "SPACE.H"
#include <cmath>
#include "FArrayBox.H"
#include "REAL.H"

#include "NamespaceHeader.H"

/// Interpolates on faces of a rectangular domain.
class InterpOnFace
/** Interpolates on faces of a rectangular domain.
*/
{

public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor leaves InterpOnFace undefined.
  */
  InterpOnFace();

  ///
  /**
     Full constructor.
     Calls full define() function with identical arguments.
  */
  InterpOnFace(/// NODE-centered face
               const Box&   a_fineFaceBox,
               /// normal direction of the face
               int          a_idir,
               /// coarsening ratio
               int          a_coarsening,
               /// width of border around face
               int          a_interpBorder);

  ///
  /** Destructor.
   */
  virtual ~InterpOnFace();

  ///
  /**
     Full define function.
  */
  virtual void define(/// NODE-centered face
                      const Box&   a_fineFaceBox,
                      /// normal direction of the face
                      int          a_idir,
                      /// coarsening ratio
                      int          a_coarsening,
                      /// width of border around face
                      int          a_interpBorder);

  /*@}*/

  /**
     \name Access functions
  */
  /*@{*/

  ///
  /** Returns <tt>true</tt> if this object was created with the defining
      constructor or if define() has been called.
  */
  bool isDefined() const;


  /*@}*/

  /**
     \name Parameter-setting functions
  */
  /*@{*/

  ///
  /** Set verbosity.  Default is 0.
   */
  void setVerbose(int a_verbose);

  /*@}*/

  /**
     \name Data modification functions
  */
  /*@{*/

  ///
  /** Interpolates from coarse to fine data on the face.
  */
  void interpolate(/// fine data
                   FArrayBox&         a_fine,
                   /// coarse data
                   const FArrayBox&   a_coarse);

  /*@}*/

protected:

  /** has this InterpOnFace been defined yet?
   */
  bool m_isDefined;

  /** NODE-centered fine face
   */
  Box m_fineFaceBox;

  /** NODE-centered coarse face
   */
  Box m_coarseFaceBox;

  /** NODE-centered coarse face refined by m_coarsening
   */
  Box m_crseFaceBoxBasicRefined;

  /** whether m_fineFaceBox - m_fineCoarseShift is coarsenable by m_coarsening
   */
  bool m_coarsenedExact;

  /** amount by which to shift a_fine after interpolation
   */
  IntVect m_offLoExt;

  /** verbosity; default 0
   */
  int m_verbose;

  /** coarsening ratio
   */
  int m_coarsening;

  /** normal direction of the face
   */
  int m_idir;

  /** parallel directions of the face
   */
  int m_ipar1, m_ipar2;

  /** width of border around face
   */
  int m_interpBorder;

  /** lower corner of m_fineFaceBox
   */
  IntVect m_fineCoarseShift;

  /** weights for interpolating to nodes where first face coordinate
      is not divisible by coarsening ratio, but second is divisible
   */
  FArrayBox m_weights1d1;

  /** weights for interpolating to nodes where first face coordinate
      is divisible by coarsening ratio, but second is not divisible
   */
  FArrayBox m_weights1d2;

  // new weights

  Box m_cbx1, m_cbx2, m_icbx1, m_cbx1e, m_cbx1n, m_icbx2;

  Box m_ifbx1, m_ifbx2, m_cfbx2e, m_cfbx2n;

  FArrayBox m_weights1, m_weights2;

  FArrayBox m_work;

private:

  void operator = (const InterpOnFace&)
  {
  }

  InterpOnFace(const InterpOnFace&)
  {
  }
};

#include "NamespaceFooter.H"

#endif
