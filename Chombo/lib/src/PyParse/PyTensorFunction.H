#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PYTENSORFUNCTION_H_
#define _PYTENSORFUNCTION_H_

// Forward declaration for PyObject.
struct _object;

#include "TensorFunction.H"
#include "NamespaceHeader.H"

//! \class PyTensorFunction
//! This class allows one to construct vector functions that are evaluated by
//! the Python interpreter.
class PyTensorFunction: public TensorFunction
{

  public:

  //! Returns true if the given Python object can be interpreted as a
  //! vector function, false otherwise.
  static bool isValid(_object* a_pyObject);

  //! Construct a vector function from the given Python object.
  //! \param a_pyObject The Python object used to construct the vector function.
  explicit PyTensorFunction(_object* a_pyObject);

  //! Destructor.
  ~PyTensorFunction();

  // Overridden methods.
  RealTensor operator()(const RealVect& a_x, Real a_t) const;

  private:

  // The underlying Python objects.
  _object *m_func, *m_derivs;
};

#ifndef PyTensorFunction_Check
#define PyTensorFunction_Check(x) PyTensorFunction::isValid(x)
#endif

#include "NamespaceFooter.H"

#endif
