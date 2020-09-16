#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MBVECTORSTENCILELEMENT_H_
#define _MBVECTORSTENCILELEMENT_H_

#include "IntVect.H"
#include "Tuple.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
/// Class to describe a single element of a multi-block interpolation stencil.
class MBVectorStencilElement
{
public:

  /// null constructor leaves object in undefined state.
  MBVectorStencilElement();

  ///
  /**
     Constructor.
  */
  MBVectorStencilElement(const IntVect&  a_cell,
                         int             a_block,
                         Tuple<Real, SpaceDim*SpaceDim> a_weight);

  // Set weight to zero.
  MBVectorStencilElement(const IntVect&  a_cell,
                         int             a_block);

  ///
  /*
    Same as constructor.
   */
  virtual void define(const IntVect&  a_cell,
                      int             a_block,
                      Tuple<Real, SpaceDim*SpaceDim> a_weight);

  // Set weight to zero.
  virtual void define(const IntVect&  a_cell,
                      int             a_block);

  bool operator==(const MBVectorStencilElement& a_input) const;

  inline
  const IntVect& cell() const;

  inline
  int block() const;

  inline
  Tuple<Real, SpaceDim*SpaceDim> weight() const;

  inline
  void setWeight(Tuple<Real, SpaceDim*SpaceDim> a_weight);

protected:

  /// index of this stencil cell
  IntVect m_cell;

  /// block containing this stencil cell
  int m_block;

  /// weight of this cell in the stencil
  Tuple<Real, SpaceDim*SpaceDim> m_weight;
};

inline
const IntVect&
MBVectorStencilElement::cell() const
{
  return m_cell;
}

inline
int
MBVectorStencilElement::block() const
{
  return m_block;
}

inline
Tuple<Real, SpaceDim*SpaceDim>
MBVectorStencilElement::weight() const
{
  return m_weight;
}

inline
void
MBVectorStencilElement::setWeight(Tuple<Real, SpaceDim*SpaceDim> a_weight)
{
  m_weight = a_weight;
}

#include "NamespaceFooter.H"

#endif // include guard
