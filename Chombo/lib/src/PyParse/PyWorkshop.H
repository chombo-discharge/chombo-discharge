#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PYWORKSHOP_H__
#define _PYWORKSHOP_H__

#ifdef CH_USE_EB

#include "Python.h"
#include "BaseIF.H"
#include "RefCountedPtr.H"

// This handles the initialization of the Workshop module.
PyMODINIT_FUNC initWorkshop();

extern PyTypeObject PyImplicitFunction_Type;
#define PyImplicitFunction_Check(obj) (obj->ob_type == &PyImplicitFunction_Type)

// Creates a PyImplicitFunction with the given name and IF.
extern PyObject* PyImplicitFunction_FromIF(const std::string& a_name,
                                           BaseIF* a_IF);

// Provide access to the underlying C++ implicit function.
extern RefCountedPtr<BaseIF> PyImplicitFunction_IF(PyObject* a_object);

#endif

#endif
