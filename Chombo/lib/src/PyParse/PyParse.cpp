#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PyParse.H"
#include "PyBox.H"
#include "PyFunctions.H"
#ifdef CH_USE_EB
#include "PyWorkshop.H"
#include "PyIF.H"
#endif
#include "Python.h"
#include <fstream>
#include "PyScalarFunction.H"
#include "ConstantScalarFunction.H"
#include "PyVectorFunction.H"
#include "ConstantVectorFunction.H"
#include "PyTensorFunction.H"
#include "ConstantTensorFunction.H"
#include "CH_assert.H"
#include "MayDay.H"

using namespace std;

struct ParsedData
{
  map<string, PyObject*> objects;
};

// Static thingies.
static PyParse::ErrorHandlerFunc s_errorHandler = NULL;
static ParsedData* s_data = NULL;

//-----------------------------------------------------------------------
// This function implements error handling.
void
PyParse::
error(const string& a_message)
{
  char str[1024];
  if (Py_IsInitialized() && PyErr_Occurred())
  {
    // The Python interpreter has something to say.
    PyObject *type, *value, *traceback;
    PyErr_Fetch(&type, &value, &traceback);
    snprintf(str, 1024, "PyParse: %s\n%s\n", a_message.c_str(), PyString_AsString(value));
    Py_XDECREF(type);
    Py_XDECREF(value);
    Py_XDECREF(traceback);
  }
  else
    snprintf(str, 1024, "PyParse: %s", a_message.c_str()); // Just use the given message.

  // If the error handler hasn't been set, use MayDay::error.
  if (s_errorHandler == NULL)
    MayDay::Error(str);
  else
    (*s_errorHandler)(str);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// This function does the heavy lifting of translating Python data to C++.
void
PyParse::
translate(PyObject* a_localVariables)
{
  // Clear any existing data.
  if (s_data != NULL)
    delete s_data;
  s_data = new ParsedData;

  // Step through the Python dictionary.
  Py_ssize_t pos = 0;
  PyObject *key, *value;
  while (PyDict_Next(a_localVariables, &pos, &key, &value))
  {
    // All items in sequences must be nonzero in length and have items
    // of identical type, period.
    if (PySequence_Check(value))
    {
      int badIndex = -1;
      Py_ssize_t len = PySequence_Length(value);
      if (len == 0)
      {
        char msg[1024];
        snprintf(msg, 1024, "Found an empty sequence '%s'.", PyString_AsString(key));
        PyParse::error(msg);
      }

      PyObject* seq = PySequence_Fast(value, "");
      PyTypeObject* pytype = NULL;
      for (int i = 0; i < len; ++i)
      {
        PyObject *item = PySequence_Fast_GET_ITEM(seq, i);
        if (pytype == NULL)
          pytype = item->ob_type;
        else if (pytype != item->ob_type)
        {
          // Okay, okay--ints can be floats.
          if (!((pytype == &PyFloat_Type) && (PyInt_Check(item))) &&
              !((pytype == &PyInt_Type) && (PyFloat_Check(item))))
            badIndex = i;
        }
      }
      Py_DECREF(seq);
      if (badIndex != -1)
      {
        char msg[1024];
        snprintf(msg, 1024, "Item %d in sequence '%s' is not a %s.",
            badIndex, PyString_AsString(key), pytype->tp_name);
        PyParse::error(msg);
      }
    }

    // Now jot down the object.
    string name = PyString_AsString(key);
    s_data->objects[name] = value;
  }

  // We're done!
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
bool
PyParse::
parsed()
{
  return (s_data != NULL);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
PyParse::
setErrorHandler(PyParse::ErrorHandlerFunc a_handler)
{
  CH_assert(a_handler != NULL);
  s_errorHandler = a_handler;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
PyParse::
parse(const string& a_file)
{
  vector<string> blank;
  parse(a_file, blank);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
PyParse::
parse(const string& a_file,
      const vector<string>& a_commands,
      const vector<string>& a_paramNames)
{
  // If we've not initialized Python yet, do so here.
  if (!Py_IsInitialized())
  {
    // Is PYTHONPATH set? If so, annihilate it.
    if (getenv("PYTHONPATH") != NULL)
      unsetenv("PYTHONPATH");

    // Set PYTHONHOME. This is necessary for the prefix to be
    // set properly on some systems.
    setenv("PYTHONHOME", CHOMBO_HOME "/src/PyParse/python", 1);

    // Initialize the Box type.
    if (PyImport_AppendInittab((char*)"BoxTools", initBoxTools) < 0)
      error("Could not initialize the Python interpreter.");

    // Tell Python about the Functions module.
    if (PyImport_AppendInittab((char*)"Functions", initFunctions) < 0)
      error("Could not initialize the Python interpreter.");
#if CH_USE_EB
    // Tell Python about the Workshop module.
    if (PyImport_AppendInittab((char*)"Workshop", initWorkshop) < 0)
      error("Could not initialize the Python interpreter.");
#endif

    // Initialize the Python interpreter.
    Py_Initialize();
    if (PyErr_Occurred())
      error("Could not initialize the Python interpreter.");

  }

  // Assemble a string containing the Python code in the file,
  // prepended by some commands to import other modules.
  ifstream file(a_file.c_str());
  if (!file)
  {
    char msg[1024];
    snprintf(msg, 1024, "Could not open input file '%s'.", a_file.c_str());
    error(msg);
  }
  string commands("from math import *\n");
  commands += string("from BoxTools import *\n");
  commands += string("from Functions import *\n");
#if CH_USE_EB
  commands += string("from Workshop import *\n");
#endif
  while (!file.eof())
  {
    char line[1024];
    file.getline(line, 1024);
    commands += line;
    commands += '\n';
  }
  file.close();
  // Add any additional commands.
  for (int i = 0; i < a_commands.size(); ++i)
    commands += a_commands[i] + string("\n");
  int status = PyRun_SimpleString(commands.c_str());
  if (status < 0)
    error("Error parsing input file.");

  // Now that we've parsed things, translate the data to C++.
  PyObject* __main__ = PyImport_ImportModule("__main__");
  PyObject* locals = PyModule_GetDict(__main__);
  if (locals == NULL)
  {
    Py_DECREF(__main__);
    error("Error obtaining Python variables from interpreter.");
  }
  translate(locals);

  if (!a_paramNames.empty())
  {
    for (int i = 0; i < a_paramNames.size(); ++i)
    {
      if (s_data->objects.find(a_paramNames[i]) == s_data->objects.end())
      {
        Py_DECREF(__main__);
        error(string("parameter ") + a_paramNames[i] + string(" not found"));
      }
    }
  }

  Py_DECREF(__main__);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
PyParse::
parse(const string& a_file,
      const vector<string>& a_paramNames)
{
  vector<string> commands;
  parse(a_file, commands, a_paramNames);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
PyParse::
parse(int a_argc,
      char* a_argv[])
{
  vector<string> blank;
  parse(a_argc, a_argv, blank);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
PyParse::
parse(int a_argc,
      char* a_argv[],
      const vector<string>& a_paramNames)
{
  // Filename comes first, followed by additional commands.
  if (a_argc < 2)
    error("PyParse: Name of input file required.");
  string filename(a_argv[1]);
  vector<string> commands;
  for (int i = 2; i < a_argc; ++i)
    commands.push_back(string(a_argv[i]));
  parse(filename, commands, a_paramNames);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyParse::
PyParse()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyParse::
PyParse(const string& a_file)
{
  parse(a_file);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyParse::
PyParse(const string& a_file,
        const vector<string>& a_paramNames)
{
  parse(a_file, a_paramNames);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyParse::
PyParse(int a_argc,
        char* a_argv[])
{
  parse(a_argc, a_argv);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyParse::
PyParse(int a_argc,
        char* a_argv[],
        const vector<string>& a_paramNames)
{
  parse(a_argc, a_argv, a_paramNames);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyParse::
PyParse(const PyParse& a_rhs)
{
  // Nothing to do here.
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyParse::
~PyParse()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PyParse&
PyParse::
operator=(const PyParse& a_rhs)
{
  return *this;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
bool
PyParse::
contains(const string& a_name)
{
  if (!parsed())
    error("PyParse has not yet successfully parsed any\ninput parameters.");
  CH_assert(s_data != NULL);
  return (s_data->objects.find(a_name) != s_data->objects.end());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
PyParse::
extract(PyObject* a_object, string& a_value)
{
  if (!PyString_Check(a_object))
    return -1;
  a_value = string(PyString_AsString(a_object));
  return 0;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
PyParse::
extract(PyObject* a_object, int& a_value)
{
  if (!PyInt_Check(a_object))
    return -1;
  a_value = PyInt_AsLong(a_object);
  return 0;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
PyParse::
extract(PyObject* a_object, bool& a_value)
{
  if (!PyInt_Check(a_object))
    return -1;
  a_value = (a_object == Py_True);
  return 0;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
PyParse::
extract(PyObject* a_object, Real& a_value)
{
  if (!PyInt_Check(a_object) && !PyFloat_Check(a_object))
    return -1;
  a_value = PyFloat_AsDouble(a_object);
  return 0;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
PyParse::
extract(PyObject* a_object, IntVect& a_value)
{
#if CH_SPACEDIM == 1
  if (!PyInt_Check(a_object) && !PySequence_Check(a_object))
    return -1;
  if (PyInt_Check(a_object))
  {
    extract(a_object, a_value[0]);
    return 0;
  }
#else
  if (!PySequence_Check(a_object))
    return -1;
#endif
  Py_ssize_t len = PySequence_Length(a_object);
  if (len < SpaceDim)
    return -1;
  PyObject* seq = PySequence_Fast(a_object, "");
  for (Py_ssize_t i = 0; i < SpaceDim; ++i)
  {
    PyObject* item = PySequence_Fast_GET_ITEM(seq, i);
    if (extract(item, a_value[i]))
    {
      Py_DECREF(seq);
      return -1;
    }
  }
  Py_DECREF(seq);
  return 0;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
PyParse::
extract(PyObject* a_object, RealVect& a_value)
{
#if CH_SPACEDIM == 1
  if (!PyFloat_Check(a_object) && !PyInt_Check(a_object) && !PySequence_Check(a_object))
    return -1;
  if (!PySequence_Check(a_object))
  {
    extract(a_object, a_value[0]);
    return 0;
  }
#else
  if (!PySequence_Check(a_object))
    return -1;
#endif
  Py_ssize_t len = PySequence_Length(a_object);
  if (len < SpaceDim)
    return -1;
  PyObject* seq = PySequence_Fast(a_object, "");
  for (Py_ssize_t i = 0; i < SpaceDim; ++i)
  {
    PyObject* item = PySequence_Fast_GET_ITEM(seq, i);
    if (extract(item, a_value[i]))
    {
      Py_DECREF(seq);
      return -1;
    }
  }
  Py_DECREF(seq);
  return 0;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
PyParse::
extract(PyObject* a_object, RealTensor& a_value)
{
#if CH_SPACEDIM == 1
  if (!PyFloat_Check(a_object) && !PyInt_Check(a_object) && !PySequence_Check(a_object))
    return -1;
  if (!PySequence_Check(a_object))
  {
    extract(a_object, a_value(0,0));
    return 0;
  }
#else
  if (!PySequence_Check(a_object))
    return -1;
#endif
  Py_ssize_t len = PySequence_Length(a_object);
  if (len < SpaceDim*SpaceDim)
    return -1;
  PyObject* seq = PySequence_Fast(a_object, "");
  for (Py_ssize_t i = 0; i < SpaceDim; ++i)
  {
    PyObject* item = PySequence_Fast_GET_ITEM(seq, i);
    int r = i / SpaceDim;
    int c = i % SpaceDim;
    if (extract(item, a_value(r,c)))
    {
      Py_DECREF(seq);
      return -1;
    }
  }
  Py_DECREF(seq);
  return 0;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
PyParse::
extract(PyObject* a_object, Box& a_value)
{
  if (!PyBox_Check(a_object))
    return -1;
  a_value = PyBox_AsBox(a_object);
  return 0;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
PyParse::
extract(PyObject* a_object, RefCountedPtr<PyUnaryFunction>& a_value)
{
  if (!PyUnaryFunction_Check(a_object))
    return -1;
  a_value = RefCountedPtr<PyUnaryFunction>(new PyUnaryFunction(a_object));
  return 0;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
PyParse::
extract(PyObject* a_object, RefCountedPtr<PyBinaryFunction>& a_value)
{
  if (!PyBinaryFunction_Check(a_object))
    return -1;
  a_value = RefCountedPtr<PyBinaryFunction>(new PyBinaryFunction(a_object));
  return 0;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
PyParse::
extract(PyObject* a_object, RefCountedPtr<ScalarFunction>& a_value)
{
  if (!PyScalarFunction_Check(a_object))
  {
    Real val;
    if (extract(a_object, val) < 0)
      return -1;
    a_value = RefCountedPtr<ScalarFunction>(new ConstantScalarFunction(val));
    return 0;
  }
  a_value = RefCountedPtr<ScalarFunction>(new PyScalarFunction(a_object));
  return 0;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
PyParse::
extract(PyObject* a_object, RefCountedPtr<VectorFunction>& a_value)
{
  if (!PyVectorFunction_Check(a_object))
  {
    RealVect val;
    if (extract(a_object, val) < 0)
      return -1;
    a_value = RefCountedPtr<VectorFunction>(new ConstantVectorFunction(val));
    return 0;
  }
  a_value = RefCountedPtr<VectorFunction>(new PyVectorFunction(a_object));
  return 0;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
PyParse::
extract(PyObject* a_object, RefCountedPtr<TensorFunction>& a_value)
{
  if (!PyTensorFunction_Check(a_object))
  {
    RealTensor val;
    if (extract(a_object, val) < 0)
      return -1;
    a_value = RefCountedPtr<TensorFunction>(new ConstantTensorFunction(val));
    return 0;
  }
  a_value = RefCountedPtr<TensorFunction>(new PyTensorFunction(a_object));
  return 0;
}
//-----------------------------------------------------------------------

#ifdef CH_USE_EB
//-----------------------------------------------------------------------
int
PyParse::
extract(PyObject* a_object, RefCountedPtr<BaseIF>& a_value)
{
  if (!PyIF_Check(a_object))
    return -1;
  a_value = RefCountedPtr<BaseIF>(new PyIF(a_object));
  return 0;
}
//-----------------------------------------------------------------------
#endif

//-----------------------------------------------------------------------
PyObject*
PyParse::
getObject(const std::string& a_name)
{
  return s_data->objects[a_name];
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
PyParse::
getSequence(PyObject* a_sequence, vector<PyObject*>& a_value)
{
  if (!PySequence_Check(a_sequence))
    return -1;
  Py_ssize_t len = PySequence_Length(a_sequence);
  a_value.resize(len);
  for (int i = 0; i < len; ++i)
  {
    a_value[i] = PySequence_GetItem(a_sequence, i);
    Py_DECREF(a_value[i]);
  }
  return 0;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
PyParse::
finalize()
{
  // Clean up our data.
  if (s_data != NULL)
  {
    delete s_data;
    s_data = NULL;
  }

  // Stop the Python parser.
  if (Py_IsInitialized())
    Py_Finalize();
}
//-----------------------------------------------------------------------

