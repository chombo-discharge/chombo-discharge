%module CHOMBO
%{
#include "IntVect.H"
#include "Box.H"
#include "ProblemDomain.H"
#include "IntVectSet.H"
#include "DataIndex.H"
#include "LayoutIterator.H"
#include "DataIterator.H"
#include "BoxLayout.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "LDF.H"
#include "AMRIO.H"
#include "Python_Utils.H"
%}
%typemap(python,in) string {
  $target = new string(PyString_AsString($source));
}
%typemap(python,in) std::string {
  $target =  new std::string(PyString_AsString($source));
}
%typemap(python,out) std::string {
  $target = PyString_FromString(($source)->c_str());
}
%typemap(python,out) string {
  $target = PyString_FromString(($source)->c_str());
}
%module outarg

      // This tells SWIG to treat an double * argument with name 'OutValue' as
      // an output value.  We'll append the value to the current result which 
      // is guaranteed to be a List object by SWIG.

 %typemap(python,argout) double *OutValue {
       PyObject *o;
              o = PyFloat_FromDouble(*$source);
              if ((!$target) || ($target == Py_None)) {
                      $target = o;
              } else {
                      if (!PyList_Check($target)) {
                              PyObject *o2 = $target;
                              $target = PyList_New(0);
                              PyList_Append($target,o2);
                              Py_XDECREF(o2);
                      }
                      PyList_Append($target,o);
                      Py_XDECREF(o);
              }
      }
%apply double { Real }; 

%include "IntVect.i"
%include "Box.i"
%include "ProblemDomain.i"
%include "IntVectSet.i"
%include "DataIndex.i"
%include "LayoutIterator.i"
%include "DataIterator.i"
%include "BoxLayout.i"
%include "DisjointBoxLayout.i"
%include "FArrayBox.i"
%include "LDF.i"
%include "AMRIO.i"
%include "Python_Utils.i"
