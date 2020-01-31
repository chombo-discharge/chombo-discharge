#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CFARRAYBOX_H_
#define _CFARRAYBOX_H_

#include "Box.H"
#include "BaseFab.H"

#include "CH_Complex.H"
#include "NamespaceHeader.H"

///
/**
   Complex Fortran Array Boxes (generally called CFArrayBoxs) are objects constructed
   to interface with arrays in Fortran.  Useful operations can be performed
   upon FABs in C++, and they provide a convenient interface to
   Fortran when it is necessary to retreat into that language for
   doing arithmetic operations when performance matters.

   CFArrayBox is derived from BaseFab<Complex>.
   CFArrayBox adds additional useful capabilities which make sense
   for Complex types, such as I/O and L**p norms.

   The C pre-processor macro `CH_SPACEDIM' must be defined to use
   this class.  The internal precision of CFArrayBox objects is
   set by defining either `CH_USE_FLOAT' or `CH_USE_DOUBLE'

   This class does NOT provide a copy constructor or assignment operator.
*/
class CFArrayBox: public BaseFab<Complex>
{
public:

  ///
  /**
     Constructs an invalid CFArrayBox with no memory.
  */
  CFArrayBox ();

  ///
  /**
     Constructs an initial CFArrayBox with the data space allocated but not
     inititialized. a_ncomp is the number of components (variables) at each
     data point in the Box.
  */
  CFArrayBox(const Box& a_box,
            int        a_ncomp,
            Complex*      a_alias = NULL);

  ///
  /**
     Construct an aliaed CFArrayBox.  See BaseFab class for details.
  */
  CFArrayBox(const Interval& a_comps,
            CFArrayBox&      a_original)
    :
    BaseFab<Complex>(a_comps, a_original)
  {}

  virtual void define(const Box& a_box, int a_ncomp, Complex* a_alias=NULL)
  {
    BaseFab<Complex>::define(a_box, a_ncomp, a_alias);
  }
  ///
  /**
     Defines CFArrayBox with the data space allocated but not
     inititialized. a_ncomp is the number of components (variables) at each
     data point in the Box.
  */
  // virtual void define(const Box& a_box,
  //                      int        a_ncomp,
  //                      Complex*      a_alias = NULL)
  //  {
  //   BaseFab<Complex>::resize(a_box, a_ncomp, a_alias);
  // }

  ///
  /**
     The (virtual) destructor.
  */
  virtual ~CFArrayBox ();

  CFArrayBox(CFArrayBox&& a_in)=default;

  CFArrayBox& operator=(CFArrayBox&& a_in)=default;
  ///
  /**
     Constructs an 'aliased' BaseFab of the requested interval of the
     argument BaseFab.  This BaseFab does not allocate any memory, but
     sets its data pointer into the memory pointed to by the argument
     BaseFab.  It is the users responsiblity to ensure this aliased
     BaseFab is not used after the original BaseFab has deleted its data ptr
     (resize, define(..) called, or destruction, etc.).

     This aliased BaseFab will also generate side effects (modifying the values
     of data in one will modify the other's data).

     This aliased BaseFab will have a_comps.size() components, starting at zero.
  */
  virtual void define(const Interval& a_comps,
                      CFArrayBox&      a_original)
  {
    BaseFab<Complex>::define(a_comps, a_original);
  }

  ///
  /**
   * This is here only to make the Intel compiler stop warning about
   * partial override.
   */
  virtual void define(const Interval& a_comps,
              BaseFab<Complex>&      a_original)
  {
    BaseFab<Complex>::define(a_comps, a_original);
  }

  ///
  /**
     Constructs CFArrayBox by reading it from istream.
  */
  explicit CFArrayBox(std::istream& a_is);

  ///
  /**
     Returns the Lp-norm of this FAB using components
     (a_comp : a_comp+a_numcomp-1) and within the a_subbox.
     a_p < 0 -> ERROR
     a_p = 0 -> infinity norm (max absolute value)
     a_p = 1 -> sum of ABS(FAB)
     a_p > 1 -> Lp-norm
  */
  virtual Real norm(const Box& a_subbox,
                    int        a_p = 2,
                    int        a_comp = 0,
                    int        a_numcomp = 1) const;

  ///
  /**
     Returns the Lp-norm of this FAB using components
     (a_comp : a_comp+a_numcomp-1).
     a_p < 0  -> ERROR
     a_p = 0  -> infinity norm (max absolute value)
     a_p = 1  -> sum of ABS(FAB)
     a_p > 1  -> Lp-norm
  */
  virtual Real norm(int a_p = 2,
                    int a_comp = 0,
                    int a_numcomp = 1) const;

  ///
  /**
     Returns sum of pow(fab[i,c],p): i in a_subbox, a_comp <= c <
     a_comp+a_numcomp, a_p >= 2 only
  */
  virtual Complex sumPow(const Box& a_subbox,
                      int        a_p = 2,
                      int        a_comp = 0,
                      int        a_numcomp = 1) const;

  /// Return the dot product of this CFArrayBox with another
  /**
     Return the dot product of this CFArrayBox and "a_fab2" over their common
     box and all components.
  */
  Complex dotProduct(const CFArrayBox& a_fab2) const;

  /// Return the dot product of this CFArrayBox with another
  /**
     Return the dot product of this CFArrayBox and "a_fab2" over the
     a_box box and all components.
  */
  Complex dotProduct(const CFArrayBox& a_fab2, const Box& a_box) const;


  ///
  /**
     Modifies this CFArrayBox by replacing each value with its absolute value.
  */
  void abs();

  ///
  /**
     Modifies this CFArrayBox by replacing each value with its absolute value,
     for components (a_comp : a_comp+a_numcomp-1).
  */
  void abs(int a_comp,
           int a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by replacing eahc value with its absolute value,
     for components (a_comp : a_comp+a_numcomp-1) and within the a_subbox.
  */
  void abs (const Box& a_subbox,
            int        a_comp = 0,
            int        a_numcomp = 1);

  ///
  /**
     Returns sum of given component of CFArrayBox.
  */
  Complex sum(int a_comp,
           int a_numcomp = 1) const;

  ///
  /**
     Returns sum of component of this CFArrayBox in given a_subbox.
  */
  Complex sum(const Box& a_subbox,
           int        a_comp,
           int        a_numcomp = 1) const;

  ///
  /**
     Modifies this CFArrayBox by replacing each value x with a_r/x.
  */
  CFArrayBox& invert(Complex a_r);

  ///
  /**
     Modifies this CFArrayBox by replacing each value x with a_r/x.  For
     given range of components.
  */
  CFArrayBox& invert(Complex a_r,
                    int  a_comp,
                    int  a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by replacing each value x with a_r/x.  For
     given range of components and within given a_subbox.
  */
  CFArrayBox& invert(Complex       a_r,
                    const Box& a_subbox,
                    int        a_comp = 0,
                    int        a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by replacing each value with its additive
     inverse.  For given range of components and within given a_subbox.
  */
  CFArrayBox& negate(const Box& a_subbox,
                    int        a_comp = 0,
                    int        a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by replacing each value with its additive
     inverse.  For given range of components.
  */
  CFArrayBox& negate(int a_comp,
                    int a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by replacing each value with its additive
     inverse.
  */
  CFArrayBox& negate();

  ///
  /**
     Modifies this CFArrayBox by adding the scalar Complex a_r to all values.  For
     given range of components and within given a_subbox.
  */
  CFArrayBox& plus(Complex       a_r,
                  const Box& a_subbox,
                  int        a_comp = 0,
                  int        a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by adding the scalar Complex a_r to all values.  For
     given range of components.
  */
  CFArrayBox& plus(Complex a_r,
                  int  a_comp,
                  int  a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by adding the scalar Complex a_r to all values.
  */
  CFArrayBox& operator += (Complex a_r);

  ///
  /**
     Modifies this CFArrayBox by pointwise addition of the values of the
     argument CFArrayBox.  You might come to grief if the domains of the
     CFArrayBoxes don't match.
  */
  CFArrayBox& operator += (const CFArrayBox& a_x);

  ///
  /**
     Modifies this CFArrayBox by adding the scalar Complex a_r to all values.
  */
  CFArrayBox& plus(Complex a_r);


  ///
  /**
     Modifies this CFArrayBox by pointwise addition of the values of the
     argument CFArrayBox.  You might come to grief if the domains of the
     CFArrayBoxes don't match.  The same as += operator.

  */
  CFArrayBox& plus(const CFArrayBox& a_x);

  ///
  /**
     Modifies this CFArrayBox by pointwise scaled addition of the
     argument CFArrayBox (a[i] <- a[i] + a_scale * a_src[i]).  Uses domain
     of the intersection of these two CFArrayBoxes.
  */
  CFArrayBox& plus(const CFArrayBox& a_src,
                  const Complex&      a_scale);

  ///
  /**
     Modifies this CFArrayBox by pointwise scaled addition of the
     argument CFArrayBox (a[i] <- a[i] + a_scale * a_src[i]).  Uses domain
     of the intersection of these two CFArrayBoxes.
  */
  CFArrayBox& plus(const CFArrayBox& a_src,
                  const Complex&      a_scale,
                  int              a_srccomp,
                  int              a_destcomp,
                  int              a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by pointwise addition of values in the argument
     CFArrayBox.  Adds src's components (a_srccomp : a_srccomp+a_numcomp-1)
     to this CFArrayBox's components (a_destcomp : a_destcomp+a_numcomp-1)
     where the domains of the two CFArrayBoxes intersect.
  */
  CFArrayBox& plus(const CFArrayBox& a_src,
                  int              a_srccomp,
                  int              a_destcomp,
                  int              a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by pointwise addition of values in the argument
     CFArrayBox.  Adds src's components (a_srccomp : a_srccomp+a_numcomp-1)
     to this CFArrayBox's components (a_destcomp : a_destcomp+numcomp-1)
     where the domain of this CFArrayBox intersects the a_subbox.  NOTE:
     a_subbox must be contained in this FAB.
  */
  CFArrayBox& plus(const CFArrayBox& a_src,
                  const Box&       a_subbox,
                  int              a_srccomp,
                  int              a_destcomp,
                  int              a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by pointwise addition of values in the argument
     CFArrayBox.  Adds src's components (a_srccomp : a_srccomp+a_numcomp-1)
     in the Box a_srcbox to this CFArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) in the Box a_destbox.  Corresponding locations
     within the two CFArrayBoxes are indexed relative to a_srcbox and a_destbox,
     and will in general not be the same.  The a_srcbox and a_destbox must be
     same size.  The results are UNDEFINED if the a_src and dest CFArrayBoxes
     are the same and the a_srcbox and a_destbox overlap.

  */
  CFArrayBox& plus(const CFArrayBox& a_src,
                  const Box&       a_srcbox,
                  const Box&       a_destbox,
                  int              a_srccomp,
                  int              a_destcomp,
                  int              a_numcomp = 1);

  ///
  CFArrayBox& plus(const CFArrayBox& a_src,
                  const Box&       a_srcbox,
                  const Box&       a_destbox,
                  const Complex&      a_scale,
                  int              a_srccomp,
                  int              a_destcomp,
                  int              a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by subtracting the scalar Complex a_r to all values.
     Note: use plus(-a_r) for more general operations.
  */
  CFArrayBox& operator -= (Complex a_r);

  ///
  /**
     Modifies this CFArrayBox by pointwise subtraction of the values of the
     argument CFArrayBox.  You might come to grief if the domains of the
     CFArrayBoxes don't match.
  */
  CFArrayBox& operator -= (const CFArrayBox& a_x);

  ///
  /**
     Modifies this CFArrayBox by pointwise subtraction of the values of the
     argument CFArrayBox.  You might come to grief if the domains of the
     CFArrayBoxes don't match.  The same as -= operator.
  */
  CFArrayBox& minus(const CFArrayBox& a_x);

  /**
     Modifies this CFArrayBox by pointwise subtraction of values in the
     argument CFArrayBox.  Subtracts a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) from this CFArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) where the domains of the two CFArrayBoxes
     intersect.
  */
  CFArrayBox& minus(const CFArrayBox& a_src,
                   int              a_srccomp,
                   int              a_destcomp,
                   int              a_numcomp = 1);

  /**
     Modifies this CFArrayBox by pointwise subtraction of values in the
     argument CFArrayBox.  Subtracts a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) from this CFArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) where the domain of this CFArrayBox intersects
     the a_subbox.  NOTE: a_subbox must be contained in this FAB.
  */
  CFArrayBox& minus(const CFArrayBox& a_src,
                   const Box&       a_subbox,
                   int              a_srccomp,
                   int              a_destcomp,
                   int              a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by pointwise subtraction of values in the
     argument CFArrayBox.  Subtracts a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) in the Box a_srcbox from this CFArrayBox's
     components (a_destcomp : a_destcomp+a_numcomp-1) in the Box a_destbox.
     Corresponding locations within the two CFArrayBoxes are indexed relative
     to a_srcbox and a_destbox, and will in general not be the same.  The
     a_srcbox and a_destbox must be same size.  The results are UNDEFINED
     if the a_src and dest CFArrayBoxes are the same and the a_srcbox and
     a_destbox overlap.
  */
  CFArrayBox& minus(const CFArrayBox& a_src,
                   const Box&       a_srcbox,
                   const Box&       a_destbox,
                   int              a_srccomp,
                   int              a_destcomp,
                   int              a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by multiplying all values by the scalar Complex a_r.
  */
  CFArrayBox& operator *= (Complex a_r);

  ///
  /**
     Modifies this CFArrayBox by multiplying all values by the scalar Complex a_r.
  */
  CFArrayBox& mult(Complex a_r);

  ///
  /**
     Modifies this CFArrayBox by multiplying all values by the scalar
     Complex a_r.  For given range of components.
  */
  CFArrayBox& mult(Complex a_r,
                  int  a_comp,
                  int  a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by multiplying all values by the scalar
     Complex a_r.  For given range of components and within given a_subbox.
  */
  CFArrayBox& mult(Complex       a_r,
                  const Box& a_subbox,
                  int        a_comp = 0,
                  int        a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by pointwise multiplication of the values by the
     argument CFArrayBox.  You might come to grief if the domains of the
     CFArrayBoxes don't match.
  */
  CFArrayBox& operator *= (const CFArrayBox& a_x);

  ///
  /**
     Modifies this CFArrayBox by pointwise multiplication by the values in the
     argument CFArrayBox.  You might come to grief if the domains of the
     CFArrayBoxes don't match.  The same as the *= operator.
  */
  CFArrayBox& mult(const CFArrayBox& a_x);

  ///
  /**
     Modifies this CFArrayBox by pointwise multiplication by values in the
     argument CFArrayBox.  Multiplies a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) by this CFArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) where the domains of the two CFArrayBoxes
     intersect.
  */
  CFArrayBox& mult(const CFArrayBox& a_src,
                  int              a_srccomp,
                  int              a_destcomp,
                  int              a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by pointwise multiplication by values in the
     argument CFArrayBox.  Multiplies a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) by this CFArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) where the domain of this CFArrayBox intersects
     the a_subbox.  NOTE: a_subbox must be contained in this FAB.
  */
  CFArrayBox& mult(const CFArrayBox& a_src,
                  const Box&       a_subbox,
                  int              a_srccomp,
                  int              a_destcomp,
                  int              a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by pointwise multiplication by values in the
     argument CFArrayBox.  Multiplies a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) in the Box a_srcbox by this CFArrayBox's
     components (a_destcomp : a_destcomp+a_numcomp-1) in the Box a_destbox.
     Corresponding locations within the two CFArrayBoxes are indexed relative
     to a_srcbox and a_destbox, and will in general not be the same.  The
     a_srcbox and a_destbox must be same size.  The results are UNDEFINED if
     the a_src and dest CFArrayBoxes are the same and the a_srcbox and a_destbox
     overlap.
  */
  CFArrayBox& mult(const CFArrayBox& a_src,
                  const Box&       a_srcbox,
                  const Box&       a_destbox,
                  int              a_srccomp,
                  int              a_destcomp,
                  int              a_numcomp = 1);

  ///
  /**
     Modify this CFArrayBox by multiplying each element by the Real value in a_r
  */
  CFArrayBox& mult(const BaseFab<Real>& a_r);

///
  /**
     Modify this CFArrayBox by multiplying each element by the Real value in a_r
     in the region defined by a_box
  */
  CFArrayBox& mult(const BaseFab<Real>& a_r,
                   const Box&           a_box);

  ///
  /**
     Modifies this CFArrayBox by dividing all values by the scalar Complex a_r.
  */
  CFArrayBox& operator /= (Complex a_r);

  ///
  /**
     Modifies this CFArrayBox by dividing all values by the scalar Complex a_r.
  */
  CFArrayBox& divide(Complex a_r);

  ///
  /**
     Modifies this CFArrayBox by dividing all values by the scalar Complex a_r.
     For given range of components.
  */
  CFArrayBox& divide(Complex a_r,
                    int  a_comp,
                    int  a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by dividing all values by the scalar Complex
     a_r.  For given range of components and within given a_subbox.
  */
  CFArrayBox& divide(Complex       a_r,
                    const Box& a_subbox,
                    int        a_comp = 0,
                    int        a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by pointwise division of the values by the
     argument CFArrayBox.  You might come to grief if the domains of the
     CFArrayBoxes don't match.
  */
  CFArrayBox& operator /= (const CFArrayBox& a_x);

  ///
  /**
     Modifies this CFArrayBox by pointwise division by the values in the
     argument CFArrayBox.  You might come to grief if the domains of the
     CFArrayBoxes don't match.  The same as the /= operator.
  */
  CFArrayBox& divide(const CFArrayBox& a_x);

  ///
  /**
     Modifies this CFArrayBox by pointwise division by values in the argument
     CFArrayBox.  Divides this CFArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) by a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) where the domains of the two CFArrayBoxes intersect.
  */
  CFArrayBox& divide(const CFArrayBox& a_src,
                    int              a_srccomp,
                    int              a_destcomp,
                    int              a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by pointwise division by values in the argument
     CFArrayBox.  Divides this CFArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) by a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) where the domain of this CFArrayBox intersects
     the a_subbox.  NOTE: a_subbox must be contained in this FAB.
  */
  CFArrayBox& divide(const CFArrayBox& a_src,
                    const Box&       a_subbox,
                    int              a_srccomp,
                    int              a_destcomp,
                    int              a_numcomp = 1);

  ///
  /**
     Modifies this CFArrayBox by pointwise division by values in the argument
     CFArrayBox.  Divides this CFArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) in the Box a_destbox by a_src's components
     (a_srccomp : a_srccomp+a_numcomp-1) in the Box a_srcbox.  Corresponding
     locations within the two CFArrayBoxes are indexed relative to a_srcbox and
     a_destbox, and will in general not be the same.  The a_srcbox and
     a_destbox must be same size.  The results are UNDEFINED if the a_src and
     dest CFArrayBoxes are the same and the a_srcbox and a_destbox overlap.
  */
  CFArrayBox& divide(const CFArrayBox& a_src,
                    const Box&       a_srcbox,
                    const Box&       a_destbox,
                    int              a_srccomp,
                    int              a_destcomp,
                    int              a_numcomp = 1);

  ///
  Complex get(const IntVect& a_iv,
           int            a_comp) const
  {
    return this->operator()(a_iv, a_comp);
  }

  ///
  void set(const IntVect& a_iv,
           int            a_comp,
           Complex           a_val)
  {
    this->operator()(a_iv, a_comp) = a_val;
  }

  //! Computes a_A * a_X + a_B * a_Y, placing the result in this CFArrayBox.
  CFArrayBox& axby(const CFArrayBox& a_X, const CFArrayBox& a_Y,
                  Complex a_A, Complex a_B);

  //! Computes a_A * a_X + a_B * a_Y, placing the result in this CFArrayBox.
  //! This version performs this operation only for the given component
  //! in each CFArrayBox.
  CFArrayBox& axby(const CFArrayBox& a_X, const CFArrayBox& a_Y,
                  Complex a_A, Complex a_B,
                  int a_destComp, int a_xComp, int a_yComp);

protected:
  virtual void performCopy(const BaseFab<Complex>& a_src,
                           const Box&           a_srcbox,
                           int                  a_srccomp,
                           const Box&           a_destbox,
                           int                  a_destcomp,
                           int                  a_numcomp);


private:
  //
  // These are disallowed.
  //
  CFArrayBox (const CFArrayBox&);
  CFArrayBox& operator = (const CFArrayBox&);
};

#include "NamespaceFooter.H"
#endif
