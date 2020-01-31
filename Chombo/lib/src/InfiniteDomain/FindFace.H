#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FINDFACE_H_
#define _FINDFACE_H_

#include "BoxLayout.H"
#include "Tuple.H"

#include "NamespaceHeader.H"

class FindFace : public BaseTransform
{
public:

  /**
     \name Constructor, destructor, and defines.
  */
  /*@{*/

  ///
  /** Constructor.
   */
  FindFace(const Box& a_baseBox,
           const Tuple<Box, SpaceDim>& a_faces);

  /** destructor
   */
  virtual ~FindFace();

  /** Required for BaseTransform.
      Returns the Box corresponding to the given face.
  */
  virtual Box operator() (const Box& a_inputBox);

protected:

  Box m_baseBox;

  Tuple<Box, SpaceDim> m_faces;
};

#include "NamespaceFooter.H"

#endif // include guard
