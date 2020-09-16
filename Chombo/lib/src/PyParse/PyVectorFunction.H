#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PYVECTORFUNCTION_H_
#define _PYVECTORFUNCTION_H_

// Forward declaration for PyObject.
struct _object;

#include "VectorFunction.H"
#include "NamespaceHeader.H"

//! \class PyVectorFunction
//! This class allows one to construct vector functions that are evaluated by
//! the Python interpreter.
class PyVectorFunction: public VectorFunction
{

  public:

  //! Returns true if the given Python object can be interpreted as a
  //! vector function, false otherwise.
  static bool isValid(_object* a_pyObject);

  //! Returns true if the given Python objects can be interpreted as a
  //! scalar function, false otherwise. The first object evaluates the
  //! function, and the second evaluates partial derivatives.
  static bool isValid(_object* a_pyFunction,
                      _object* a_pyDerivs);

  //! Construct a vector function from the given Python object.
  //! \param a_pyObject The Python object used to construct the vector function.
  explicit PyVectorFunction(_object* a_pyObject);

  //! Construct a vector function from given Python objects. This
  //! constructor provides a method for evaluating partial derivatives
  //! of the function.
  //! \param a_pyFunction The Python object used to evaluate the vector function.
  //! \param a_pyDerivs The Python object used to evaluate the partial derivatives.
  PyVectorFunction(_object* a_pyFunction,
                   _object* a_pyDerivs);

  //! Destructor.
  ~PyVectorFunction();

  // Overridden methods.
  RealVect operator()(const RealVect& a_x, Real a_t) const;
  RealVect derivative(const IntVect& a_order, const RealVect& a_x, Real a_t) const;
  bool hasDerivative(const IntVect& a_order) const;

  private:

  // The underlying Python objects.
  _object *m_func, *m_derivs;
};

#ifndef PyVectorFunction_Check
#define PyVectorFunction_Check(x) PyVectorFunction::isValid(x)
#endif

#include "NamespaceFooter.H"

#endif
