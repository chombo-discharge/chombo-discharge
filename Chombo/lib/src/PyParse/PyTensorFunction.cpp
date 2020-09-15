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
#include "PyTensorFunction.H"
#include "CH_assert.H"
#include <iostream>
#include "NamespaceHeader.H"

using namespace std;

//-----------------------------------------------------------------------
bool
PyTensorFunction::
isValid(PyObject* a_pyObject)
{
  // If we are given a tensor, we can interpret the function as constant
  // in both space and time.
  if (PySequence_Check(a_pyObject) && (PySequence_Length(a_pyObject) == SpaceDim*SpaceDim))
  {
    bool isTensor = true;
    for (int d = 0; d < SpaceDim*SpaceDim; ++d)
    {
      PyObject* item = PySequence_GetItem(a_pyObject, d);
      if (!PyFloat_Check(item) && !PyInt_Check(item))
        isTensor = false;
      Py_DECREF(item);
    }
    return isTensor;
  }

  // Otherwise, if it can't be called, it is an abomination.
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

  // The result should be interpretable as a tensor.
  bool isTensor = (PySequence_Check(result) &&
                   (PySequence_Length(result) == SpaceDim*SpaceDim));
  if (isTensor)
  {
    for (int d = 0; d < SpaceDim*SpaceDim; ++d)
    {
      PyObject* item = PySequence_GetItem(result, d);
      if (!PyFloat_Check(item) && !PyInt_Check(item))
        isTensor = false;
      Py_DECREF(item);
    }
  }
  Py_DECREF(result);
  return isTensor;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyTensorFunction::
PyTensorFunction(PyObject* a_pyObject):
  TensorFunction(false, false),
  m_func(a_pyObject),
  m_derivs(NULL)
{
  CH_assert(a_pyObject != NULL);
  CH_assert(isValid(a_pyObject));

  // If we are given a tensor, we can interpret the function as constant
  // in both space and time.
  if (PySequence_Check(a_pyObject) && (PySequence_Length(a_pyObject) == SpaceDim*SpaceDim))
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
    // The 2-argument form doesn't work, so the function is constant in time.
    PyErr_Clear();
    m_isConstant = true;
  }
  Py_XDECREF(result);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyTensorFunction::
~PyTensorFunction()
{
  Py_XDECREF(m_func);
  Py_XDECREF(m_derivs);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
RealTensor
PyTensorFunction::
operator()(const RealVect& a_x, Real a_t) const
{
  // If we're constant and homogeneous, interpret the function
  // as a tensor.
  RealTensor t;
  if (m_isConstant && m_isHomogeneous)
  {
    for (int i = 0; i < SpaceDim; ++i)
    {
      for (int j = 0; j < SpaceDim; ++j)
      {
        PyObject* item = PyTuple_GetItem(m_func, SpaceDim*j+i);
        t(i, j) = (PyFloat_Check(item)) ? static_cast<double>(PyFloat_AsDouble(item))
                                        : static_cast<double>(PyInt_AsLong(item));
      }
    }
    return t;
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
    MayDay::Error("Evaluation of Python-supplied vector function failed.");
  }

  // Extract our vector.
  CH_assert(PySequence_Check(result) && (PySequence_Length(result) == SpaceDim));
  for (int i = 0; i < SpaceDim; ++i)
  {
    for (int j = 0; j < SpaceDim; ++j)
    {
      PyObject* item = PySequence_GetItem(result, SpaceDim*j+i);
      t(i,j) = (PyFloat_Check(item)) ? static_cast<double>(PyFloat_AsDouble(item))
                                     : static_cast<double>(PyInt_AsLong(item));
      Py_DECREF(item);
    }
  }
  Py_DECREF(result);
  return t;
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
