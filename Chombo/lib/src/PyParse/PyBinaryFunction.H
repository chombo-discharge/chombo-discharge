#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PYBINARYFUNCTION_H_
#define _PYBINARYFUNCTION_H_

// Forward declaration for PyObject.
struct _object;

#include "REAL.H"

//! \class PyBinaryFunction
//! This class allows one to construct binary functions that are evaluated by
//! the Python interpreter. A binary function is a function \f$f\f$ that
//! takes two (numeric) arguments \f$x\f$ and \f$y\f$, returning a single
//! real-valued number \f$f(x,y)\f$.
class PyBinaryFunction
{

  public:

  //! Returns true if the given Python object can be interpreted as a
  //! Binary function, false otherwise. Such an object is a binary function if
  //! it is able to accept two Reals or ints and turn them into another float or int.
  //! \param a_pyObject The Python object in question.
  static bool isValid(_object* a_pyObject);

  //! Construct a Binary function from the given Python object.
  //! \param a_pyObject The Python object used to construct the Binary function.
  explicit PyBinaryFunction(_object* a_pyObject);

  //! Destructor.
  ~PyBinaryFunction();

  //! Evaluation operator. This calls the corresponding Python function, returning
  //! the result as a Real.
  //! \param a_x The first argument to the function.
  //! \param a_x The second argument to the function.
  Real operator()(Real a_x, Real a_y) const;

  private:

  // The underlying Python object.
  _object* m_obj;

  // Disallowed stuff.
  PyBinaryFunction();
  PyBinaryFunction(const PyBinaryFunction&);
  PyBinaryFunction& operator=(const PyBinaryFunction&);
};

#ifndef PyBinaryFunction_Check
#define PyBinaryFunction_Check(x) PyBinaryFunction::isValid(x)
#endif

#endif
