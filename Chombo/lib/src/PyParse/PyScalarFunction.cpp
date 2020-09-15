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
#include "PyScalarFunction.H"
#include "CH_assert.H"
#include "NamespaceHeader.H"
using namespace std;

//-----------------------------------------------------------------------
bool
PyScalarFunction::
isValid(PyObject* a_pyObject)
{
  // If we are given a number, we can interpret the function as constant
  // in both space and time.
  if (PyFloat_Check(a_pyObject) || PyInt_Check(a_pyObject))
    return true;

  // If the object is a 2-tuple with valid contents, it's also okay.
  if (PyTuple_Check(a_pyObject) &&
      (PyTuple_Size(a_pyObject) == 2) &&
      PyScalarFunction::isValid(PyTuple_GetItem(a_pyObject, 0),
                                PyTuple_GetItem(a_pyObject, 1)))
    return true;

  // Try to call the object with a single argument and then watch what
  // happens.

  // First of all, if it can't be called, it is an abomination.
  if (!PyCallable_Check(a_pyObject))
    return false;

  // Try to call it at time 0 at the origin.
  PyObject* x = PyTuple_New(SpaceDim);
  for (int d = 0; d < SpaceDim; ++d)
  {
    PyTuple_SetItem(x, d, PyFloat_FromDouble(1e-15));
  }
  double t = 0.0;
  PyObject* result = PyObject_CallFunction(a_pyObject, (char*)"Od", x, t);
  if (result == NULL)
  {
    // The 2-argument form doesn't work. Try the single-argument form.
    PyErr_Clear();
    result = PyObject_CallFunctionObjArgs(a_pyObject, x, NULL);
    if (result == NULL)
    {
      // Nope. No good.
      PyErr_Clear();
      Py_DECREF(x);
      return false;
    }
  }
  Py_DECREF(x);

  // The result should be interpretable as a number.
  bool isNumber = (PyFloat_Check(result) || PyInt_Check(result));
  Py_DECREF(result);
  return isNumber;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
bool
PyScalarFunction::
isValid(PyObject* a_pyFunction,
        PyObject* a_pyDerivs)
{
  // First of all, if they can't be called, they are abominations.
  if (!PyCallable_Check(a_pyFunction) || !PyCallable_Check(a_pyDerivs))
    return false;

  // Try to call the function at time 0 at the origin.
  PyObject* order = PyTuple_New(SpaceDim);
  PyObject* x = PyTuple_New(SpaceDim);
  for (int d = 0; d < SpaceDim; ++d)
  {
    PyTuple_SetItem(order, d, PyInt_FromLong(1));
    PyTuple_SetItem(x, d, PyFloat_FromDouble(1e-15));
  }
  double t = 0.0;
  PyObject* derivResult;
  PyObject* result = PyObject_CallFunction(a_pyFunction, (char*)"Od", x, t);
  if (result == NULL)
  {
    // The 2-argument form doesn't work. Try the single-argument form.
    PyErr_Clear();
    result = PyObject_CallFunctionObjArgs(a_pyFunction, x, NULL);
    if (result == NULL)
    {
      // Nope. No good.
      PyErr_Clear();
      Py_DECREF(order);
      Py_DECREF(x);
      return false;
    }
    // Okay, the single-argument form of the function works, which means
    // that it is constant in time. Therefore, the 2-argument form of
    // the derivative should work.
    derivResult = PyObject_CallFunction(a_pyDerivs, (char*)"OO", order, x);
    if (derivResult == NULL)
    {
      // Well, we almost made it.
      PyErr_Clear();
      Py_DECREF(order);
      Py_DECREF(x);
      return false;
    }
  }
  else
  {
    // The 2-argument form of the function works, which means
    // that it is time-dependent. Therefore, the 3-argument form of
    // the derivative should work.
    derivResult = PyObject_CallFunction(a_pyDerivs, (char*)"OOd", order, x, t);
    if (derivResult == NULL)
    {
      // No luck.
      PyErr_Clear();
      Py_DECREF(order);
      Py_DECREF(x);
      return false;
    }
  }
  Py_DECREF(order);
  Py_DECREF(x);

  // The results should be interpretable as numbers.
  bool funcIsNumber = (PyFloat_Check(result) || PyInt_Check(result));
  bool derivIsNumber = (PyFloat_Check(derivResult) || PyInt_Check(derivResult));
  Py_DECREF(result);
  Py_DECREF(derivResult);
  return (funcIsNumber && derivIsNumber);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyScalarFunction::
PyScalarFunction(PyObject* a_pyObject):
  ScalarFunction(false, false),
  m_func(a_pyObject),
  m_derivs(NULL)
{
  CH_assert(a_pyObject != NULL);
  CH_assert(isValid(a_pyObject));

  // If the object is a 2-tuple with valid contents, crack it open.
  if (PyTuple_Check(a_pyObject))
  {
    m_func = PyTuple_GetItem(a_pyObject, 0);
    m_derivs = PyTuple_GetItem(a_pyObject, 1);
  }
  Py_XINCREF(m_func);
  Py_XINCREF(m_derivs);

  // If the object is a number, the function is constant.
  if (PyFloat_Check(a_pyObject) || PyInt_Check(a_pyObject))
  {
    m_isConstant = m_isHomogeneous = true;
    return;
  }

  // Figure out whether the function takes one or two arguments.
  PyObject* x = PyTuple_New(SpaceDim);
  for (int d = 0; d < SpaceDim; ++d)
  {
    PyTuple_SetItem(x, d, PyFloat_FromDouble(1e-15));
  }
  double t = 0.0;
  PyObject* result = PyObject_CallFunction(a_pyObject, (char*)"Od", x, t);
  Py_DECREF(x);
  if (result == NULL)
  {
    // The 2-argument form doesn't work, so the function is constant.
    PyErr_Clear();
    m_isConstant = true;
  }
  Py_XDECREF(result);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyScalarFunction::
PyScalarFunction(PyObject* a_pyFunction,
                 PyObject* a_pyDerivs):
  ScalarFunction(false, false),
  m_func(a_pyFunction),
  m_derivs(a_pyDerivs)
{
  CH_assert(a_pyFunction != NULL);
  CH_assert(a_pyDerivs != NULL);
  CH_assert(isValid(a_pyFunction, a_pyDerivs));
  Py_INCREF(a_pyFunction);
  Py_INCREF(a_pyDerivs);

  // Figure out whether the function takes one or two arguments.
  PyObject* x = PyTuple_New(SpaceDim);
  for (int d = 0; d < SpaceDim; ++d)
  {
    PyTuple_SetItem(x, d, PyFloat_FromDouble(1e-15));
  }
  double t = 0.0;
  PyObject* result = PyObject_CallFunction(a_pyFunction, (char*)"Od", x, t);
  Py_DECREF(x);
  if (result == NULL)
  {
    // The 2-argument form doesn't work, so the function is constant.
    PyErr_Clear();
    m_isConstant = true;
  }
  Py_XDECREF(result);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyScalarFunction::
~PyScalarFunction()
{
  Py_XDECREF(m_func);
  Py_XDECREF(m_derivs);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
PyScalarFunction::
operator()(const RealVect& a_x, Real a_t) const
{
  // If we're constant and homogeneous, interpret the function
  // as a number.
  if (m_isConstant && m_isHomogeneous)
  {
    Real val = (PyFloat_Check(m_func)) ? static_cast<double>(PyFloat_AsDouble(m_func))
                                       : static_cast<double>(PyInt_AsLong(m_func));
    return val;
  }

  // Construct a d-tuple for a_x.
  PyObject* x = PyTuple_New(SpaceDim);
  for (int d = 0; d < SpaceDim; ++d)
  {
    PyTuple_SetItem(x, d, PyFloat_FromDouble(a_x[d]));
  }

  // Call the function with our arguments.
  PyObject* result;
  if (isConstant())
  {
    result = PyObject_CallFunctionObjArgs(m_func, x, NULL);
  }
  else
  {
    result = PyObject_CallFunction(m_func, (char*)"Od", x,
                                   static_cast<double>(a_t));
  }
  Py_DECREF(x);
  if (PyErr_Occurred())
  {
    PyErr_Print();
    MayDay::Error("Evaluation of Python-supplied scalar function failed.");
  }
  Real val = (PyFloat_Check(result)) ? static_cast<double>(PyFloat_AsDouble(result))
                                     : static_cast<double>(PyInt_AsLong(result));
  Py_DECREF(result);
  return val;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
PyScalarFunction::
derivative(const IntVect& a_order,
           const RealVect& a_x,
           Real a_t) const
{
  // If the order is zero, return the function's value.
  if (a_order == IntVect::Zero)
    return (*this)(a_x, a_t);

  // If the function is homogeneous, return 0 for all derivatives.
  if (m_isHomogeneous)
    return 0.0;

  // If we don't have a mechanism for evaluating derivatives,
  // we are unable to continue.
  if (m_derivs == NULL)
    MayDay::Error("This scalar function does not define derivatives");

  // Construct d-tuples for a_order and a_x.
  PyObject* order = PyTuple_New(SpaceDim);
  PyObject* x = PyTuple_New(SpaceDim);
  for (int d = 0; d < SpaceDim; ++d)
  {
    PyTuple_SetItem(order, d, PyInt_FromLong(a_order[d]));
    PyTuple_SetItem(x, d, PyFloat_FromDouble(a_x[d]));
  }

  // Call the function with our arguments.
  PyObject* result;
  if (isConstant())
  {
    result = PyObject_CallFunction(m_derivs, (char*)"OO", order, x);
  }
  else
  {
    result = PyObject_CallFunction(m_derivs, (char*)"OOd", order, x,
                                   static_cast<double>(a_t));
  }
  Py_DECREF(order);
  Py_DECREF(x);
  if (PyErr_Occurred())
  {
    PyErr_Print();
    MayDay::Error("Evaluation of Python-supplied scalar function derivative failed.");
  }
  Real val = (PyFloat_Check(result)) ? static_cast<double>(PyFloat_AsDouble(result))
                                     : static_cast<double>(PyInt_AsLong(result));
  Py_DECREF(result);
  return val;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
bool
PyScalarFunction::
hasDerivative(const IntVect& a_order) const
{
  // We always have the zeroth order derivative. :-P
  if ((a_order == IntVect::Zero) || m_isHomogeneous)
    return true;

  // If we don't have a mechanism for evaluating derivatives,
  // we are unable to continue.
  if (m_derivs == NULL)
    return false;

  // Construct d-tuples for a_order and a_x.
  PyObject* order = PyTuple_New(SpaceDim);
  PyObject* x = PyTuple_New(SpaceDim);
  for (int d = 0; d < SpaceDim; ++d)
  {
    PyTuple_SetItem(order, d, PyInt_FromLong(a_order[d]));
    PyTuple_SetItem(x, d, PyFloat_FromDouble(1e-15));
  }

  // Call the function with our arguments.
  PyObject* result;
  if (isConstant())
  {
    result = PyObject_CallFunction(m_derivs, (char*)"OO", order, x);
  }
  else
  {
    result = PyObject_CallFunction(m_derivs, (char*)"OOd", order, x,
                                   static_cast<double>(0.0));
  }
  Py_DECREF(order);
  Py_DECREF(x);
  if (PyErr_Occurred())
  {
    PyErr_Clear();
    return false;
  }
  Py_DECREF(result);
  return true;
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
