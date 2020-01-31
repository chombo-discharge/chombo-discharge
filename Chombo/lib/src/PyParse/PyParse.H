#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// By Jeffrey Johnson (resident Python zealot)

#ifndef _PYPARSE_H_
#define _PYPARSE_H_

#include <string>
#include <vector>
#include <map>
#include "Vector.H"
#include "RefCountedPtr.H"
#include "REAL.H"
#include "IntVect.H"
#include "RealVect.H"
#include "Box.H"
#include "PyUnaryFunction.H"
#include "PyBinaryFunction.H"
#include "ScalarFunction.H"
#include "VectorFunction.H"
#include "TensorFunction.H"
#ifdef CH_USE_EB
#include "BaseIF.H"
#endif

// Forward declaration for PyObject.
struct _object;

//! \class PyParse
//! This class parses parameters from a Python file. It is similar to
//! ParmParse, but since it uses the Python interpreter it is at once
//! simpler and more capable.
class PyParse
{
  public:

  //! \typedef ErrorHandlerFunc
  //! An ErrorHandlerFunc is a void function that accepts an error
  //! message as input.
  typedef void(*ErrorHandlerFunc)(const char*);

  //! Returns true if PyParse has parsed a file yet, false otherwise.
  static bool parsed();

  //! Sets the error handler for the parser. By default, this is set to
  //! MayDay::error.
  //! \param a_handler A function that is called when an error occurs.
  static void setErrorHandler(ErrorHandlerFunc a_handler);

  //! Parse the given file, gathering all parameters therein and storing
  //! them in C++ types.
  //! \param a_file The name of a file containing Python code.
  static void parse(const std::string& a_file);

  //! Parse the given string, interpreting it as a Python statement.
  //! \param a_command The name of a file containing Python code.
  static void parseCommand(const std::string& a_command);

  //! Parse the given file, gathering only the parameters whose names appear
  //! in the given vector. If \a a_paramNames is empty, all parameters found
  //! are parsed.
  //! \param a_file The name of a file containing Python code.
  //! \param a_paramNames A vector containing names of the parameters to
  //!                     retrieve.
  static void parse(const std::string& a_file,
                    const std::vector<std::string>& a_paramNames);

  //! Parse the given file and any additional commands.
  //! \param a_file The name of a file containing Python code.
  //! \param a_commands A vector containing additional commands.
  //! \param a_paramNames A vector containing names of the parameters to
  //!                     retrieve.
  static void parse(const std::string& a_file,
                    const std::vector<std::string>& a_commands,
                    const std::vector<std::string>& a_paramNames);

  //! Use argc and argv to find the given file and parse it. This
  //! command also parses any command line assignments it encounters
  //! and overwrites the corresponding variables in the file.
  //! \param a_argc The number of command line arguments.
  //! \param a_argv The array of command line arguments.
  static void parse(int a_argc, char* a_argv[]);

  //! Use argc and argv to find the given file and parse it, gathering only the
  //! parameters whose names appear in the given vector. If \a a_paramNames is
  //! empty, all parameters found are parsed.
  //! \param a_argc The number of command line arguments.
  //! \param a_argv The array of command line arguments.
  //! \param a_paramNames A vector containing names of the parameters to
  //!                     retrieve.
  static void parse(int a_argc, char* a_argv[],
                    const std::vector<std::string>& a_paramNames);

  //! Default constructor. This provides access to the parser and can be
  //! queried if parsing has already been done.
  PyParse();

  //! Parse the given file, gathering all parameters therein and storing them
  //! in C++ types. This constructor creates a PyParse object, which is just
  //! a proxy to the parsed data.
  //! \param a_file The name of a file containing Python code.
  explicit PyParse(const std::string& a_file);

  //! Parse the given file, gathering only the parameters whose names appear
  //! in the given vector. This constructor creates a PyParse object, which is just
  //! a proxy to the parsed data. If \a a_paramNames is empty, all parameters
  //! found are parsed.
  //! \param a_file The name of a file containing Python code.
  //! \param a_paramNames A vector containing names of the parameters to
  //!                     retrieve.
  PyParse(const std::string& a_file,
          const std::vector<std::string>& a_paramNames);

  //! Use argc and argv to find the given file and parse it.  This constructor
  //! creates a PyParse object, which is just a proxy to the parsed data.
  //! \param a_argc The number of command line arguments.
  //! \param a_argv The array of command line arguments.
  PyParse(int a_argc, char* a_argv[]);

  //! Use argc and argv to find the given file and parse it, gathering only the
  //! parameters whose names appear in the given vector. This constructor
  //! creates a PyParse object, which is just a proxy to the parsed data.
  //! If \a a_paramNames is empty, all parameters found are parsed.
  //! \param a_argc The number of command line arguments.
  //! \param a_argv The array of command line arguments.
  //! \param a_paramNames A vector containing names of the parameters to
  //!                     retrieve.
  PyParse(int a_argc, char* a_argv[],
          const std::vector<std::string>& a_paramNames);

  //! Copy constructor.
  PyParse(const PyParse& a_rhs);

  //! Destructor.
  ~PyParse();

  //! Assignment operator.
  PyParse& operator=(const PyParse& a_rhs);

  //! Returns true if PyParse has parsed any parameter with the given name.
  //! \param a_name The name of the parameter.
  static bool contains(const std::string& a_name);

  //! Retrieves the parameter with the given name and type. If this parameter is
  //! not found, the error handler is called.
  //! \param a_name The name of the parameter.
  //! \param a_value The parameter's value will be stored here.
  template <typename T>
  static void get(const std::string& a_name, T& a_value)
  {
    _object* obj = getObject(a_name);
    if (obj == NULL)
    {
      error(a_name + std::string(" was not found!"));
    }

    // Try extracting the object itself.
    if (extract(obj, a_value) < 0)
    {
      // If that didn't work, it may be a sequence. If so, try
      // extracting the first object and using that.
      std::vector<_object*> seq;
      if (getSequence(obj, seq) < 0)
        error(a_name + " could not be properly interpreted.");
      CH_assert(!seq.empty());
      if (extract(seq[0], a_value) < 0)
        error(a_name + " could not be properly interpreted.");
    }
  }

  //! Query function for backward compatibility. This function
  //! checks to see whether the given quantity exists within the
  //! parser, and retrieves its value if it does. If it does
  //! not exist, this function has no effect.
  template <typename T>
  static bool query(const std::string& a_name, T& a_value)
  {
    bool result = contains(a_name);
    if (result)
      get(a_name, a_value);
    return result;
  }

  //! Call this to break down the Python parser.
  static void finalize();

  private:

  //! Access to parsed objects and their items.
  static _object* getObject(const std::string& a_name);
  static int getSequence(_object* a_sequence, std::vector<_object*>& a_value);

  //! Extraction methods. To add a new C++/Chombo data type to PyParse, simply
  //! implement a corresponding extract method below. Vectors of these types are
  //! automatically supported by template voodoo.
  static int extract(_object* a_object, std::string& a_value);
  static int extract(_object* a_object, int& a_value);
  static int extract(_object* a_object, bool& a_value);
  static int extract(_object* a_object, Real& a_value);
  static int extract(_object* a_object, IntVect& a_value);
  static int extract(_object* a_object, RealVect& a_value);
  static int extract(_object* a_object, RealTensor& a_value);
  static int extract(_object* a_object, Box& a_value);
  static int extract(_object* a_object, RefCountedPtr<PyUnaryFunction>& a_value);
  static int extract(_object* a_object, RefCountedPtr<PyBinaryFunction>& a_value);
  static int extract(_object* a_object, RefCountedPtr<ScalarFunction>& a_value);
  static int extract(_object* a_object, RefCountedPtr<VectorFunction>& a_value);
  static int extract(_object* a_object, RefCountedPtr<TensorFunction>& a_value);
#ifdef CH_USE_EB
  static int extract(_object* a_object, RefCountedPtr<BaseIF>& a_value);
#endif

  //! Generic vector extraction.
  template <typename T>
  static int extract(_object* a_object, std::vector<T>& a_value)
  {
    std::vector<_object*> seq;
    if (getSequence(a_object, seq) < 0)
    {
      // If a_object isn't a sequence, try interpreting it as a T.
      a_value.resize(1);
      if (extract(a_object, a_value[0]) < 0)
        return -1;
      return 0;
    }
    a_value.resize(seq.size());
    for (int i = 0; i < a_value.size(); ++i)
    {
      _object* obj = seq[i];
      if (extract(obj, a_value[i]) < 0)
        return -1;
    }
    return 0;
  }

  //! Support for Chombo Vector class.
  template <typename T>
  static int extract(_object* a_object, Vector<T>& a_value)
  {
    return extract(a_object, a_value.stdVector());
  }

  //! Error reporting.
  static void error(const std::string& a_message);

  //! Translation.
  static void translate(_object* a_localVariables);

};

#endif
