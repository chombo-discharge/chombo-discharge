#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifdef CH_USE_EB

#include "Python.h" // This must be first to side-step a Python 2.7 bug that prevents C++ code from building.
#include <string>
#include "RefCountedPtr.H"
#include "IndexTM.H"
#include "UnionIF.H"
#include "IntersectionIF.H"
#include "ComplementIF.H"
#include "PlaneIF.H"
#include "SphereIF.H"
#include "LatheIF.H"
#include "EllipsoidIF.H"
#include "TorusIF.H"
#include "PyWorkshop.H"
using namespace std;

//-----------------------------------------------------------------------
static int
getRealVect(RealVect& x, PyObject* pyx, const string& identifier)
{
  CH_assert(PyTuple_Check(pyx));
  if (PyTuple_Size(pyx) != SpaceDim)
  {
    PyErr_Format(PyExc_TypeError, "%s must be a %d-tuple.", identifier.c_str(),
                 SpaceDim);
    return -1;
  }

  for (int d = 0; d < SpaceDim; ++d)
  {
    PyObject* xd = PyTuple_GetItem(pyx, d);
    if (!PyInt_Check(xd) && !PyFloat_Check(xd))
    {
      PyErr_Format(PyExc_TypeError,
        "Item %d in %s is not a number!", d, identifier.c_str());
      return -1;
    }
    x[d] = (PyFloat_Check(xd)) ? static_cast<Real>(PyFloat_AsDouble(xd))
                               : static_cast<Real>(PyInt_AsLong(xd));
  }
  return 0;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static int
getIndexTM(IndexTM<Real, GLOBALDIM>& x, PyObject* pyx, const string& identifier)
{
  CH_assert(PyTuple_Check(pyx));
  if (PyTuple_Size(pyx) != GLOBALDIM)
  {
    PyErr_Format(PyExc_TypeError, "%s must be a %d-tuple.", identifier.c_str(),
                 GLOBALDIM);
    return -1;
  }

  for (int d = 0; d < GLOBALDIM; ++d)
  {
    PyObject* xd = PyTuple_GetItem(pyx, d);
    if (!PyInt_Check(xd) && !PyFloat_Check(xd))
    {
      PyErr_Format(PyExc_TypeError,
        "Item %d in %s is not a number!", d, identifier.c_str());
      return -1;
    }
    x[d] = (PyFloat_Check(xd)) ? static_cast<Real>(PyFloat_AsDouble(xd))
                               : static_cast<Real>(PyInt_AsLong(xd));
  }
  return 0;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static const char*
IndexTM_repr(const IndexTM<Real, GLOBALDIM>& x)
{
  static char repr[1024];
  if (GLOBALDIM == 2)
    snprintf(repr, 1024, "(%g, %g)", x[0], x[1]);
  else if (GLOBALDIM == 3)
    snprintf(repr, 1024, "(%g, %g, %g)", x[0], x[1], x[2]);
  else if (GLOBALDIM == 4)
    snprintf(repr, 1024, "(%g, %g, %g, %g)", x[0], x[1], x[2], x[3]);
  else if (GLOBALDIM == 5)
    snprintf(repr, 1024, "(%g, %g, %g, %g, %g)", x[0], x[1], x[2], x[3], x[4]);
  else
    snprintf(repr, 1024, "(%g, %g, %g, %g, %g, %g)", x[0], x[1], x[2], x[3], x[4], x[5]);
  return (const char*)repr;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyDoc_STRVAR(PyImplicitFunction_doc,
  "ImplicitFunction -- A python interface to implicit functions.\n\n");
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
typedef struct
  {
    PyObject_HEAD
    string* name;
    RefCountedPtr<BaseIF>* func;
  } PyImplicitFunction;
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyObject *
PyImplicitFunction_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyImplicitFunction* func = (PyImplicitFunction*)type->tp_alloc(type, 0);
  func->name = new string();
  func->func = new RefCountedPtr<BaseIF>();
  return (PyObject*)func;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static void
PyImplicitFunction_dealloc(PyImplicitFunction* self)
{
  delete self->name;
  delete self->func;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyObject*
PyImplicitFunction_str(PyImplicitFunction* self)
{
  return PyString_FromString(self->name->c_str());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyGetSetDef PyImplicitFunction_getsetters[] =
{
  // sentinel
  {
    NULL
  }
};
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyMethodDef PyImplicitFunction_methods[] =
{
#if 0
  {"current", (PyCFunction)PyImplicitFunction_I, METH_O,
     PyDoc_STR("I = currentSource.current(t)\n"
               "Computes the total current at time t.");
  },
#endif
  // sentinel
  {
    NULL,
    NULL
  }
};
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyObject*
PyImplicitFunction_union(PyImplicitFunction* self,
                         PyObject* rhs)
{
  if (!PyImplicitFunction_Check(rhs))
  {
    PyErr_SetString(PyExc_TypeError,
      "Set operations are only allowed between two ImplicitFunctions.");
    return NULL;
  }
  PyImplicitFunction* rIF = (PyImplicitFunction*)rhs;
  string unionName = "Union(" + *self->name + ", " + *rIF->name + ")";
  BaseIF* unionFunc = new UnionIF(**self->func, **rIF->func);
  return PyImplicitFunction_FromIF(unionName, unionFunc);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyObject*
PyImplicitFunction_iunion(PyImplicitFunction* self,
                          PyObject* rhs)
{
  if (!PyImplicitFunction_Check(rhs))
  {
    PyErr_SetString(PyExc_TypeError,
      "Set operations are only allowed between two ImplicitFunctions.");
    return NULL;
  }
  PyImplicitFunction* rIF = (PyImplicitFunction*)rhs;
  string unionName = "Union(" + *self->name + ", " + *rIF->name + ")";
  BaseIF* unionFunc = new UnionIF(**self->func, **rIF->func);
  *(self->func) = RefCountedPtr<BaseIF>(unionFunc);
  Py_INCREF((PyObject*)self);
  return (PyObject*)self;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyObject*
PyImplicitFunction_intersection(PyImplicitFunction* self,
                                PyObject* rhs)
{
  if (!PyImplicitFunction_Check(rhs))
  {
    PyErr_SetString(PyExc_TypeError,
      "Set operations are only allowed between two ImplicitFunctions.");
    return NULL;
  }
  PyImplicitFunction* rIF = (PyImplicitFunction*)rhs;
  string intName = "Intersection(" + *self->name + ", " + *rIF->name + ")";
  BaseIF* intFunc = new IntersectionIF(**self->func, **rIF->func);
  return PyImplicitFunction_FromIF(intName, intFunc);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyObject*
PyImplicitFunction_complement(PyImplicitFunction* self)
{
  string compName = "~(" + *self->name + ")";
  BaseIF* compFunc = new ComplementIF(**self->func);
  return PyImplicitFunction_FromIF(compName, compFunc);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyObject*
PyImplicitFunction_difference(PyImplicitFunction* self,
                              PyObject* rhs)
{
  if (!PyImplicitFunction_Check(rhs))
  {
    PyErr_SetString(PyExc_TypeError,
      "Set operations are only allowed between two ImplicitFunctions.");
    return NULL;
  }

  // We construct the difference by intersecting this IF with the complement of rhs.
  PyImplicitFunction* rIF = (PyImplicitFunction*)rhs;
  string diffName = "Difference(" + *self->name + ", " + *rIF->name + ")";
  BaseIF* comprIF = new ComplementIF(**rIF->func);
  BaseIF* diffFunc = new IntersectionIF(**self->func, *comprIF);
  delete comprIF;
  return PyImplicitFunction_FromIF(diffName, diffFunc);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyObject*
PyImplicitFunction_idifference(PyImplicitFunction* self,
                               PyObject* rhs)
{
  if (!PyImplicitFunction_Check(rhs))
  {
    PyErr_SetString(PyExc_TypeError,
      "Set operations are only allowed between two ImplicitFunctions.");
    return NULL;
  }
  // We construct the difference by intersecting this IF with the complement of rhs.
  PyImplicitFunction* rIF = (PyImplicitFunction*)rhs;
  string diffName = "Difference(" + *self->name + ", " + *rIF->name + ")";
  BaseIF* comprIF = new ComplementIF(**rIF->func);
  *(self->func) = RefCountedPtr<BaseIF>(new IntersectionIF(**self->func, *comprIF));
  delete comprIF;
  Py_INCREF((PyObject*)self);
  return (PyObject*)self;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyObject*
PyImplicitFunction_float(PyImplicitFunction* self)
{
  PyErr_SetString(PyExc_TypeError,
    "Cannot convert an ImplicitFunction to a float.");
  return NULL;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyNumberMethods PyImplicitFunction_as_number =
{
  (binaryfunc)PyImplicitFunction_union, /*nb_add*/
  (binaryfunc)PyImplicitFunction_difference, /*nb_subtract*/
  0, /*nb_multiply*/
  0, /*nb_divide*/
  0, /*nb_remainder*/
  0, /*nb_divmod*/
  0, /*nb_power*/
  (unaryfunc)PyImplicitFunction_complement, /*nb_negative*/
  0, /*nb_positive*/
  0, /*nb_absolute*/
  0, /*nb_nonzero*/
  0, /*nb_invert*/
  0, /*nb_lshift*/
  0, /*nb_rshift*/
  (binaryfunc)PyImplicitFunction_intersection, /*nb_and*/
  (binaryfunc)PyImplicitFunction_difference, /*nb_xor*/
  (binaryfunc)PyImplicitFunction_union, /*nb_or*/
  0, /*nb_coerce*/
  0, /*nb_int*/
  0, /*nb_long*/
  (unaryfunc)PyImplicitFunction_float, /*nb_float*/
  0, /* nb_oct */
  0, /* nb_hex */
  (binaryfunc)PyImplicitFunction_iunion, /* nb_inplace_add */
  (binaryfunc)PyImplicitFunction_idifference, /* nb_inplace_subtract */
  0, /* nb_inplace_multiply */
  0, /* nb_inplace_divide */
  0, /* nb_inplace_remainder */
  0, /* nb_inplace_power */
  0, /* nb_inplace_lshift */
  0, /* nb_inplace_rshift */
  (binaryfunc)PyImplicitFunction_intersection, /* nb_inplace_and */
  (binaryfunc)PyImplicitFunction_difference, /* nb_inplace_xor */
  (binaryfunc)PyImplicitFunction_union, /* nb_inplace_or */
  0, /* nb_floor_divide */
  0, /* nb_true_divide */
  0, /* nb_inplace_floor_divide */
  0 /* nb_inplace_true_divide */
};
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyObject*
PyImplicitFunction_call(PyObject* self,
                        PyObject* args,
                        PyObject* kwds)
{
  // For now, we only support the call to the function itself--no derivs.
  PyObject *pyx;
  if (!PyArg_ParseTuple(args, "O!", &PyTuple_Type, &pyx))
    return NULL;

  if (PyTuple_Size(pyx) != SpaceDim)
  {
    PyErr_Format(PyExc_TypeError, "Argument must be a %d-tuple.", SpaceDim);
    return NULL;
  }
  RealVect x;
  for (int d = 0; d < SpaceDim; ++d)
  {
    PyObject* xd = PyTuple_GetItem(pyx, d);
    if (!PyInt_Check(xd) && !PyFloat_Check(xd))
    {
      PyErr_Format(PyExc_TypeError,
        "Item %d in argument is not a number!", d);
      return NULL;
    }
    x[d] = (PyFloat_Check(xd)) ? static_cast<Real>(PyFloat_AsDouble(xd))
                               : static_cast<Real>(PyInt_AsLong(xd));
  }

  PyImplicitFunction* func = (PyImplicitFunction*)self;
  return PyFloat_FromDouble((*func->func)->value(x));
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyTypeObject PyImplicitFunction_Type =
{
  PyObject_HEAD_INIT(NULL)
  0,                       /*ob_size*/
  "Workshop.ImplicitFunction",   /*tp_name*/
  sizeof(PyImplicitFunction), /*tp_basicsize*/
  0,      /*tp_itemsize*/
  (destructor)PyImplicitFunction_dealloc, /*tp_dealloc*/
  0,      /*tp_print*/
  0,      /*tp_getattr*/
  0,      /*tp_setattr*/
  0,      /*tp_compare*/
  (reprfunc)PyImplicitFunction_str, /*tp_repr*/
  &PyImplicitFunction_as_number, /*tp_as_number*/
  0,      /*tp_as_sequence*/
  0,      /*tp_as_mapping*/
  0,      /*tp_hash*/
  (ternaryfunc)PyImplicitFunction_call, /*tp_call*/
  (reprfunc)PyImplicitFunction_str, /*tp_str*/
  0,      /*tp_getattro*/
  0,      /*tp_setattro*/
  0,      /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT, /*tp_flags*/
  PyImplicitFunction_doc, /*tp_doc*/
  0,      /*tp_traverse*/
  0,      /*tp_clear*/
  0,      /*tp_richcompare*/
  0,      /*tp_weaklistoffset*/
  0,      /*tp_iter*/
  0,      /*tp_iternext*/
  PyImplicitFunction_methods, /*tp_methods*/
  0,      /*tp_members*/
  PyImplicitFunction_getsetters, /*tp_getset*/
  0,      /*tp_base*/
  0,      /*tp_dict*/
  0,      /*tp_descr_get*/
  0,      /*tp_descr_set*/
  0,      /*tp_dictoffset*/
  0,      /*tp_init*/
  0,      /*tp_alloc*/
  PyImplicitFunction_new,      /*tp_new*/
  0,      /*tp_free*/
  0,      /*tp_is_gc*/
};
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Creates a new PyImplicitFunction from a C++ pointer to an IF.
PyObject*
PyImplicitFunction_FromIF(const string& a_name,
                          BaseIF* a_IF)
{
  PyImplicitFunction* result =
    (PyImplicitFunction*)PyImplicitFunction_new(&PyImplicitFunction_Type, 0, 0);
  *result->name = a_name;
  *result->func = RefCountedPtr<BaseIF>(a_IF);
  return (PyObject*)result;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Provides access to the underlying C++ model.
RefCountedPtr<BaseIF>
PyImplicitFunction_IF(PyObject* a_object)
{
  return *(((PyImplicitFunction*)a_object)->func);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyDoc_STRVAR(Workshop_doc,
  "Workshop - A Python module for constructing implicit functions with\n"
  "Solid Constructive Geometry.\n\n");
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyObject*
Workshop_Plane(PyObject* module, PyObject* args, PyObject* kwds)
{
  PyObject *pyn, *pyx, *pyinside = Py_True;
  static char* kwlist[] =
  {
    (char*)"normal",
    (char*)"point",
    (char*)"inside",
    0
  };
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|O!", kwlist,
      &PyTuple_Type, &pyn, &PyTuple_Type, &pyx, &PyBool_Type, &pyinside))
    return NULL;

  RealVect n;
  if (getRealVect(n, pyn, "normal") < 0);
    return NULL;
  RealVect x;
  if (getRealVect(x, pyx, "point") < 0);
    return NULL;

  char name[1024];
#if CH_SPACEDIM == 2
  snprintf(name, 1024, "Plane(n = (%g, %g), x = (%g, %g))", n[0], n[1], x[0], x[1]);
#else
  snprintf(name, 1024, "Plane(n = (%g, %g, %g), x = (%g, %g, %g))", n[0], n[1], n[2], x[0], x[1], x[2]);
#endif
  BaseIF* func = new PlaneIF(n, x, (pyinside == Py_True));
  return PyImplicitFunction_FromIF(name, func);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyObject*
Workshop_Sphere(PyObject* module, PyObject* args, PyObject* kwds)
{
  double r;
  PyObject *pyx, *pyinside = Py_True;
  static char* kwlist[] =
  {
    (char*)"radius",
    (char*)"center",
    (char*)"inside",
    0
  };
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "dO!|O!", kwlist,
      &r, &PyTuple_Type, &pyx, &PyBool_Type, &pyinside))
    return NULL;

  if (r <= 0.0)
  {
    PyErr_SetString(PyExc_ValueError, "radius must be positive.");
    return NULL;
  }

  RealVect x;
  if (getRealVect(x, pyx, "center") < 0)
    return NULL;

  char name[1024];
#if CH_SPACEDIM == 2
  snprintf(name, 1024, "Sphere(r = %g, x = (%g, %g))", r, x[0], x[1]);
#else
  snprintf(name, 1024, "Sphere(r = %g, x = (%g, %g, %g))", r, x[0], x[1], x[2]);
#endif
  BaseIF* func = new SphereIF(static_cast<Real>(r), x, (pyinside == Py_True));
  return PyImplicitFunction_FromIF(name, func);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyObject*
Workshop_Ellipsoid(PyObject* module, PyObject* args, PyObject* kwds)
{
  PyObject *pyr, *pyx, *pyinside = Py_True;
  static char* kwlist[] =
  {
    (char*)"radii",
    (char*)"center",
    (char*)"inside",
    0
  };
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|O!", kwlist,
      &PyTuple_Type, &pyr, &PyTuple_Type, &pyx, &PyBool_Type, &pyinside))
    return NULL;

  RealVect r;
  if (getRealVect(r, pyr, "center") < 0)
    return NULL;
  RealVect x;
  if (getRealVect(x, pyx, "center") < 0)
    return NULL;

  char name[1024];
#if CH_SPACEDIM == 2
  snprintf(name, 1024, "Ellipsoid(r = (%g, %g), x = (%g, %g))", r[0], r[1], x[0], x[1]);
#else
  snprintf(name, 1024, "Ellipsoid(r = (%g, %g, %g), x = (%g, %g, %g))", r[0], r[1], r[2], x[0], x[1], x[2]);
#endif
  BaseIF* func = new EllipsoidIF(r, x, (pyinside == Py_True));
  return PyImplicitFunction_FromIF(name, func);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyObject*
Workshop_Lathe(PyObject* module, PyObject* args, PyObject* kwds)
{
  PyObject *pyf1, *pyf2 = NULL, *pyx = NULL, *pyinside = Py_True;
  static char* kwlist[] =
  {
    (char*)"f1",
    (char*)"f2",
    (char*)"point",
    (char*)"inside",
    0
  };
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|O!O!O!", kwlist,
      &PyImplicitFunction_Type, &pyf1, &PyImplicitFunction_Type, &pyf2, &PyTuple_Type,
      &PyBool_Type, &pyinside))
    return NULL;

  PyImplicitFunction *f1 = (PyImplicitFunction*)pyf1,
                     *f2 = (PyImplicitFunction*)pyf2;
  char name[1024];
  BaseIF* func;
  if (pyf2 != NULL)
  {
    if (pyx == NULL)
    {
      PyErr_SetString(PyExc_ValueError,
        "point must be specified if f2 is given.");
      return NULL;
    }
    if (PyTuple_Size(pyx) != SpaceDim)
    {
      PyErr_Format(PyExc_ValueError, "Point must be a %d-tuple.", SpaceDim);
      return NULL;
    }
    RealVect x;
    if (getRealVect(x, pyx, "center") < 0)
      return NULL;

#if CH_SPACEDIM == 2
    snprintf(name, 1024, "Lathe(f1 = %s, f2 = %s, x = (%g, %g))", f1->name->c_str(),
             f2->name->c_str(), x[0], x[1]);
#else
    snprintf(name, 1024, "Lathe(f1 = %s, f2 = %s, x = (%g, %g, %g))", f1->name->c_str(),
             f2->name->c_str(), x[0], x[1], x[2]);
#endif
    func = new LatheIF(*PyImplicitFunction_IF(pyf1), *PyImplicitFunction_IF(pyf2), x,
                       (pyinside == Py_True));
  }
  else
  {
    snprintf(name, 1024, "Lathe(f1 = %s)", f1->name->c_str());
    func = new LatheIF(*PyImplicitFunction_IF(pyf1), (pyinside == Py_True));
  }
  return PyImplicitFunction_FromIF(name, func);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyObject*
Workshop_Torus(PyObject* module, PyObject* args, PyObject* kwds)
{
  double R, r;
  PyObject *pyx, *pyinside = Py_True;
  static char* kwlist[] =
  {
    (char*)"majorRadius",
    (char*)"minorRadius",
    (char*)"center",
    (char*)"inside",
    0
  };
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "ddO!|O!", kwlist,
      &R, &r, &PyTuple_Type, &pyx, &PyBool_Type, &pyinside))
    return NULL;

  if ((r <= 0.0) || (R <= 0.0))
  {
    PyErr_SetString(PyExc_ValueError, "Both major and minor radius must be positive.");
    return NULL;
  }
  if (r > 0.5*R)
  {
    PyErr_SetString(PyExc_ValueError,
      "minorRadius must be less than 1/2 the major radius.");
    return NULL;
  }

  RealVect x;
  if (getRealVect(x, pyx, "center") < 0)
    return NULL;

  char name[1024];
#if CH_SPACEDIM == 2
  snprintf(name, 1024, "Torus(R = %g, r = %g, x = (%g, %g))", R, r, x[0], x[1]);
#else
  snprintf(name, 1024, "Sphere(R = %g, r = %g, x = (%g, %g, %g))", R, r, x[0], x[1], x[2]);
#endif
  BaseIF* func = new TorusIF(static_cast<Real>(R), static_cast<Real>(r),
                             x, (pyinside == Py_True));
  return PyImplicitFunction_FromIF(name, func);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
static PyMethodDef Workshop_methods[] =
{
  {
   "Plane", (PyCFunction)Workshop_Plane, METH_VARARGS | METH_KEYWORDS,
     "Plane(normal, point, inside = True)\n"
     "Returns an ImplicitFunction representing a plane with the given\n"
     "normal vector centered at the given point. If inside is True, the\n"
     "plane's normal points to the domain."
  },
  {
   "Ellipsoid", (PyCFunction)Workshop_Ellipsoid, METH_VARARGS | METH_KEYWORDS,
     "Ellipsoid(radii, center, inside = True)\n"
     "Returns an ImplicitFunction representing an ellipsoid. radius is a D-tuple\n"
     "holding the radii of the ellipse, and center is a D-tuple holding the\n"
     "coordinates of the center of the ellipsoid. If inside is True, the\n"
     "ellipsoid's normal points to the domain."
  },
  {
   "Sphere", (PyCFunction)Workshop_Sphere, METH_VARARGS | METH_KEYWORDS,
     "Sphere(radius, center, inside = True)\n"
     "Returns an ImplicitFunction representing a sphere with the given radius\n"
     "and center. If inside is True, the sphere's normal points to the domain."
  },
  {
   "Lathe", (PyCFunction)Workshop_Lathe, METH_VARARGS | METH_KEYWORDS,
     "Lathe(f1, inside = True) or\n"
     "Lathe(f1, f2, point, inside = True)\n"
     "Returns an ImplicitFunction that uses one or two implicit functions to\n"
     "produce a generalized surface of revolution. If inside is True, the\n"
     "surface's normal points to the domain."
  },
  {
   "Torus", (PyCFunction)Workshop_Torus, METH_VARARGS | METH_KEYWORDS,
     "Torus(majorRadius, minorRadius, center, inside = True)\n"
     "Returns an ImplicitFunction representing a torus with the given radii\n"
     "and center. If inside is True, the torus's normal points to the domain."
  },
  // sentinel
  {
    NULL,
    NULL
  }
};
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
addModuleAttributes(PyObject* module)
{
//  PyObject* pyKtoeV = PyFloat_FromDouble(KtoeV);
//  PyModule_AddObject(module, "KtoeV", pyKtoeV);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyMODINIT_FUNC
initWorkshop()
{
  // Finalize the type object including setting type of the new type
  // object; doing it here is required for portability to Windows
  // without requiring C++.
  if (PyType_Ready(&PyImplicitFunction_Type) < 0)
    return;

  // Create the module and add its methods.
  PyObject* m = Py_InitModule3("Workshop", Workshop_methods, Workshop_doc);
  if (m == NULL)
    return;

  // Add our types.
  Py_INCREF(&PyImplicitFunction_Type);
  PyModule_AddObject(m, "ImplicitFunction", (PyObject*)&PyImplicitFunction_Type);

  // Add miscellaneous attributes.
  addModuleAttributes(m);
}
//-----------------------------------------------------------------------

#endif


