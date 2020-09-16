#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PYSCALARFUNCTION_H_
#define _PYSCALARFUNCTION_H_

// Forward declaration for PyObject.
struct _object;

#include "ScalarFunction.H"
#include "NamespaceHeader.H"

//! \class PyScalarFunction
//! This class allows one to construct scalar functions that are evaluated by
//! the Python interpreter.
class PyScalarFunction: public ScalarFunction
{

  public:

  //! Returns true if the given Python object can be interpreted as a
  //! scalar function, false otherwise.
  static bool isValid(_object* a_pyObject);

  //! Returns true if the given Python objects can be interpreted as a
  //! scalar function, false otherwise. The first object evaluates the
  //! function, and the second evaluates partial derivatives.
  static bool isValid(_object* a_pyFunction,
                      _object* a_pyDerivs);

  //! Construct a scalar function from the given Python object. This
  //! constructor does not define a method for evaluating derivatives of
  //! this function.
  //! \param a_pyObject The Python object used to construct the scalar function.
  explicit PyScalarFunction(_object* a_pyObject);

  //! Construct a scalar function from given Python objects. This
  //! constructor provides a method for evaluating partial derivatives
  //! of the function.
  //! \param a_pyFunction The Python object used to evaluate the scalar function.
  //! \param a_pyDerivs The Python object used to evaluate the partial derivatives.
  PyScalarFunction(_object* a_pyFunction,
                   _object* a_pyDerivs);

  //! Destructor.
  ~PyScalarFunction();

  // Overridden methods.

  Real operator()(const RealVect& a_x, Real a_t) const;
  Real derivative(const IntVect& a_order, const RealVect& a_x, Real a_t) const;
  bool hasDerivative(const IntVect& a_order) const;

  private:

  // The underlying Python objects.
  _object *m_func, *m_derivs;
};

#ifndef PyScalarFunction_Check
#define PyScalarFunction_Check(x) PyScalarFunction::isValid(x)
#endif

#include "NamespaceFooter.H"

#endif
