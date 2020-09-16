#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PYUNARYFUNCTION_H_
#define _PYUNARYFUNCTION_H_

// Forward declaration for PyObject.
struct _object;

#include "REAL.H"

//! \class PyUnaryFunction
//! This class allows one to construct unary functions that are evaluated by
//! the Python interpreter. A unary function is a function \f$f\f$ that
//! takes a single (numeric) argument \f$x\f$ and returns a single real-valued
//! number \f$f(x)\f$.
class PyUnaryFunction
{

  public:

  //! Returns true if the given Python object can be interpreted as a
  //! unary function, false otherwise. Such an object is a unary function if
  //! it is able to accept a Real or an int and turn it into another float or int.
  //! \param a_pyObject The Python object in question.
  static bool isValid(_object* a_pyObject);

  //! Construct a unary function from the given Python object.
  //! \param a_pyObject The Python object used to construct the unary function.
  explicit PyUnaryFunction(_object* a_pyObject);

  //! Destructor.
  ~PyUnaryFunction();

  //! Evaluation operator. This calls the corresponding Python function, returning
  //! the result as a Real.
  //! \param a_x The single argument to the function.
  Real operator()(Real a_x) const;

  private:

  // The underlying Python object.
  _object* m_obj;

  // Disallowed stuff.
  PyUnaryFunction();
  PyUnaryFunction(const PyUnaryFunction&);
  PyUnaryFunction& operator=(const PyUnaryFunction&);
};

#ifndef PyUnaryFunction_Check
#define PyUnaryFunction_Check(x) PyUnaryFunction::isValid(x)
#endif

#endif
