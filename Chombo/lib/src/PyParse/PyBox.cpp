#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <sstream>
#include "Python.h"
#include "Box.H"

//-----------------------------------------------------------------------
PyDoc_STRVAR(PyBox_doc,
  "Box(lo, hi) -> A Box in index space.\n\n");
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
typedef struct
{
  PyObject_HEAD
  Box* box;
}
PyBox;
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyObject *
PyBox_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyBox* box = (PyBox*)type->tp_alloc(type, 0);
  box->box = NULL;
  return (PyObject*)box;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static void
PyBox_dealloc(PyBox* self)
{
  if (self->box != NULL)
   delete self->box;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyObject*
PyBox_str(PyBox* self)
{
  std::ostringstream str;
  str << *self->box << '\0';
  return PyString_FromString(str.str().c_str());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyGetSetDef PyBox_getsetters[] =
{
  // sentinel
  {
    NULL
  }
};
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyMethodDef PyBox_methods[] =
{
  // sentinel
  {
    NULL,
    NULL
  }
};
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static int
PyBox_init(PyObject *self,
           PyObject *args,
           PyObject *kwds)
{
  PyBox* box = (PyBox*)self;

  PyObject *pylo, *pyhi;
  static const char* kwlist[] =
  {
    "lo",
    "hi",
    0
  };
  if (!PyArg_ParseTupleAndKeywords(args, kwds,
        "OO", (char**)kwlist, &pylo, &pyhi))
    return -1;

  IntVect lo, hi;
  if (PySequence_Check(pylo) && PySequence_Length(pylo) != SpaceDim)
  {
    PyErr_Format(PyExc_TypeError, "lo must be a %d-tuple.", SpaceDim);
    return -1;
  }
  if (PySequence_Check(pyhi) && PySequence_Length(pyhi) != SpaceDim)
  {
    PyErr_Format(PyExc_TypeError, "hi must be a %d-tuple.", SpaceDim);
    return -1;
  }
#if CH_SPACEDIM == 1
  if (!PyInt_Check(pylo) && !PySequence_Check(pylo))
  {
    PyErr_SetString(PyExc_TypeError, "lo must be an integer or 1-tuple.");
    return -1;
  }
  if (!PyInt_Check(pyhi) && !PySequence_Check(pyhi))
  {
    PyErr_SetString(PyExc_TypeError, "hi must be an integer or 1-tuple.");
    return -1;
  }
  if (PyInt_Check(pylo))
    lo[0] = PyInt_AsLong(pylo);
  if (PyInt_Check(pyhi))
    hi[0] = PyInt_AsLong(pyhi);
#else
  if (!PySequence_Check(pylo))
  {
    PyErr_Format(PyExc_TypeError, "lo must be a %d-tuple.", SpaceDim);
    return -1;
  }
  if (!PySequence_Check(pyhi))
  {
    PyErr_Format(PyExc_TypeError, "hi must be a %d-tuple.", SpaceDim);
    return -1;
  }
#endif
  if (PySequence_Check(pylo))
  {
    for (int d = 0; d < SpaceDim; ++d)
    {
      PyObject *item = PySequence_GetItem(pylo, d);
      if (!PyInt_Check(item))
      {
        PyErr_Format(PyExc_TypeError, "lo[%d] must be an integer.", d);
        Py_DECREF(item);
        return -1;
      }
      lo[d] = PyInt_AsLong(item);
      Py_DECREF(item);
    }
  }
  if (PySequence_Check(pyhi))
  {
    for (int d = 0; d < SpaceDim; ++d)
    {
      PyObject *item = PySequence_GetItem(pyhi, d);
      if (!PyInt_Check(item))
      {
        PyErr_Format(PyExc_TypeError, "hi[%d] must be an integer.", d);
        Py_DECREF(item);
        return -1;
      }
      hi[d] = PyInt_AsLong(item);
      Py_DECREF(item);
    }
  }
  box->box = new Box(lo, hi);
  return 0;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyTypeObject PyBox_Type =
{
  PyObject_HEAD_INIT(NULL)
  0,                       /*ob_size*/
  "Box",   /*tp_name*/
  sizeof(PyBox), /*tp_basicsize*/
  0,      /*tp_itemsize*/
  (destructor)PyBox_dealloc, /*tp_dealloc*/
  0,      /*tp_print*/
  0,      /*tp_getattr*/
  0,      /*tp_setattr*/
  0,      /*tp_compare*/
  (reprfunc)PyBox_str, /*tp_repr*/
  0, /*tp_as_number*/
  0,      /*tp_as_sequence*/
  0,      /*tp_as_mapping*/
  0,      /*tp_hash*/
  0,      /*tp_call*/
  (reprfunc)PyBox_str, /*tp_str*/
  0,      /*tp_getattro*/
  0,      /*tp_setattro*/
  0,      /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT, /*tp_flags*/
  PyBox_doc, /*tp_doc*/
  0,      /*tp_traverse*/
  0,      /*tp_clear*/
  0,      /*tp_richcompare*/
  0,      /*tp_weaklistoffset*/
  0,      /*tp_iter*/
  0,      /*tp_iternext*/
  PyBox_methods, /*tp_methods*/
  0,      /*tp_members*/
  PyBox_getsetters, /*tp_getset*/
  0,      /*tp_base*/
  0,      /*tp_dict*/
  0,      /*tp_descr_get*/
  0,      /*tp_descr_set*/
  0,      /*tp_dictoffset*/
  PyBox_init,      /*tp_init*/
  0,      /*tp_alloc*/
  PyBox_new,      /*tp_new*/
  0,      /*tp_free*/
  0,      /*tp_is_gc*/
};
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyObject*
PyBox_FromBox(const Box& a_box)
{
  PyBox* result =
    (PyBox*)PyBox_new(&PyBox_Type, 0, 0);
  result->box = new Box(a_box);
  return (PyObject*)result;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Box
PyBox_AsBox(PyObject* a_object)
{
  return *(((PyBox*)a_object)->box);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyDoc_STRVAR(BoxTools_doc,
  "BoxTools - A Python module for constructing Boxes.\n\n");
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyMethodDef BoxTools_methods[] =
{
  // sentinel
  {
    NULL,
    NULL
  }
};
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyMODINIT_FUNC
initBoxTools()
{
  // Initialize the Box type.
  if (PyType_Ready(&PyBox_Type) < 0)
    return;

  // Create the module and add its methods.
  PyObject* m = Py_InitModule3("BoxTools", BoxTools_methods, BoxTools_doc);
  if (m == NULL)
    return;

  // Add our types.
  Py_INCREF(&PyBox_Type);
  PyModule_AddObject(m, "Box", (PyObject*)&PyBox_Type);

  // Add SpaceDim to give input scripts dimensional awareness.
  PyModule_AddObject(m, "SpaceDim", PyInt_FromLong(SpaceDim));
}
//-----------------------------------------------------------------------


