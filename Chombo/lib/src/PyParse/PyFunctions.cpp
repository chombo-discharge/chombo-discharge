#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PyFunctions.H"
#include "PyScalarFunction.H"
#include "PyVectorFunction.H"
using namespace std;

//-----------------------------------------------------------------------
PyDoc_STRVAR(Functions_doc,
  "Functions - A Python module for constructing scalar and vector\n"
  "functions for problem generation.\n\n");
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyObject*
Functions_ScalarFunction(PyObject* module, PyObject* args, PyObject* kwds)
{
  PyObject *func, *derivs;
  static char* kwlist[] =
  {
    (char*)"function",
    (char*)"derivatives",
    0
  };
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO", kwlist,
      &func, &derivs))
    return NULL;

  if (!PyScalarFunction::isValid(func, derivs))
  {
    PyErr_SetString(PyExc_TypeError,
      "The given arguments do not represent a valid scalar function.");
    return NULL;
  }

  // Return a 2-tuple. PyParse will know what to do with this.
  return Py_BuildValue("(OO)", func, derivs);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyObject*
Functions_VectorFunction(PyObject* module, PyObject* args, PyObject* kwds)
{
  PyObject *func, *derivs;
  static char* kwlist[] =
  {
    (char*)"function",
    (char*)"derivatives",
    0
  };
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO", kwlist,
      &func, &derivs))
    return NULL;

  if (!PyVectorFunction::isValid(func, derivs))
  {
    PyErr_SetString(PyExc_TypeError,
      "The given arguments do not represent a valid vector function.");
    return NULL;
  }

  // Return a 2-tuple. PyParse will know what to do with this.
  return Py_BuildValue("(OO)", func, derivs);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyMethodDef Functions_methods[] =
{
  {
    "ScalarFunction", (PyCFunction)Functions_ScalarFunction, METH_VARARGS | METH_KEYWORDS,
    "ScalarFunction(function, derivatives)\n"
    "Returns a representation of a scalar function defined by a callable\n"
    "object (whose value at x[,t] is given by function(x[,t]) and having\n"
    "spatial derivatives defined by a callable object (whose\n"
    "values at x[,t] are given by derivatives(p, x[, t]), with p a multi-index\n"
    "identifying the order of the partial derivative to be computed."
  },
  {
    "VectorFunction", (PyCFunction)Functions_VectorFunction, METH_VARARGS | METH_KEYWORDS,
    "VectorFunction(function, derivatives)\n"
    "Returns a representation of a vector function defined by a callable\n"
    "object (whose value at x[,t] is given by function(x[,t]) and having\n"
    "spatial derivatives defined by a callable object (whose\n"
    "values at x[,t] are given by derivatives(p, x[, t]), with p a multi-index\n"
    "identifying the order of the partial derivative to be computed."
  },
  // sentinel
  {
    NULL,
    NULL
  }
};
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyMODINIT_FUNC
initFunctions()
{
  // Create the module and add its methods.
  PyObject* m = Py_InitModule3("Functions", Functions_methods, Functions_doc);
  if (m == NULL)
    return;
}
//-----------------------------------------------------------------------


