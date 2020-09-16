#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Python.h"
#include "PyBinaryFunction.H"
#include "CH_assert.H"

//-----------------------------------------------------------------------
bool
PyBinaryFunction::
isValid(PyObject* a_pyObject)
{
  // Try to call the object with a single argument and then watch what
  // happens.

  // First of all, if it can't be called, it is an abomination.
  if (!PyCallable_Check(a_pyObject))
    return false;

  // Try to call it with the value 1, since that's typically a warm and fuzzy
  // value. If it fails in a way that suggests it's a domain problem, we'll
  // try to accommodate.
  // FIXME: Figure out a formal way to do this.
  double value = 1.0;
  PyObject* result = PyObject_CallFunction(a_pyObject, (char*)"dd", value, value);
  if (result == NULL)
  {
    PyErr_Clear();
    return false;
  }

  // The result should be interpretable as a number.
  bool isNumber = (PyFloat_Check(result) || PyInt_Check(result));
  Py_DECREF(result);
  return isNumber;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyBinaryFunction::
PyBinaryFunction(PyObject* a_pyObject):
  m_obj(a_pyObject)
{
  CH_assert(a_pyObject != NULL);
  CH_assert(isValid(a_pyObject));
  Py_INCREF(a_pyObject); // Grab another reference.
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyBinaryFunction::
~PyBinaryFunction()
{
  Py_XDECREF(m_obj); // Release the reference.
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
PyBinaryFunction::
operator()(Real a_x, Real a_y) const
{
  // Call the function with our argument.
  PyObject* result = PyObject_CallFunction(m_obj, (char*)"dd",
                       static_cast<double>(a_x), static_cast<double>(a_y));
  CH_assert(result != NULL);
  Real val = (PyFloat_Check(result)) ? static_cast<double>(PyFloat_AsDouble(result))
                                     : static_cast<double>(PyInt_AsLong(result));
  Py_DECREF(result);
  return val;
}
//-----------------------------------------------------------------------

