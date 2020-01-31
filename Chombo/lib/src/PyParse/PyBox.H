#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PYBOX_H_
#define _PYBOX_H_

#include "Python.h"

extern PyTypeObject PyBox_Type;
#define PyBox_Check(obj) (obj->ob_type == &PyBox_Type)

// This handles the initialization of the Functions module.
PyMODINIT_FUNC initBoxTools();

// Creates a new PyBox from a C++ Box.
PyObject* PyBox_FromBox(const Box& a_box);

// Provides access to the underlying C++ Box.
Box PyBox_AsBox(PyObject* a_object);

#endif
