#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MBSTENCILELEMENT_H_
#define _MBSTENCILELEMENT_H_

#include "IntVect.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
/// Class to describe a single element of a multi-block interpolation stencil.
class MBStencilElement
{
public:

  /// null constructor leaves object in undefined state.
  MBStencilElement();

  ///
  /**
     Constructor.
  */
  MBStencilElement(const IntVect&  a_cell,
                   int             a_block,
                   Real            a_weight = 0.);

  ///
  /*
    Same as constructor.
   */
  virtual void define(const IntVect&  a_cell,
                      int             a_block,
                      Real            a_weight = 0.);

  bool operator==(const MBStencilElement& a_input) const;

  inline
  const IntVect& cell() const;

  inline
  int block() const;

  inline
  Real weight() const;

  inline
  void setWeight(Real a_weight);

protected:

  /// index of this stencil cell
  IntVect m_cell;

  /// block containing this stencil cell
  int m_block;

  /// weight of this cell in the stencil
  Real m_weight;
};

inline
const IntVect&
MBStencilElement::cell() const
{
  return m_cell;
}

inline
int
MBStencilElement::block() const
{
  return m_block;
}

inline
Real
MBStencilElement::weight() const
{
  return m_weight;
}

inline
void
MBStencilElement::setWeight(Real a_weight)
{
  m_weight = a_weight;
}

#include "NamespaceFooter.H"

#endif // include guard
