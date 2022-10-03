.. _Chap:CodeStandard:

Code standard
=============

When submitting new code to ``chombo-discharge``, the following guidelines below show be followed.

C++ standard
------------

We are currently at ``C++14``. 

Namespace
---------

All code in ``chombo-discharge`` is embedded in a namespace ``ChomboDischarge``.
Embedding into a namespace is done by including header file :file:`CD_NamespaceHeader.H` that contain the necessary definitions.
This is done by including after any other file includes.
In addition, files must include :file:`CD_NamespaceFooter.H` at the end. 

File names
----------

Each file should contain only one class definition, and the file name must be name of the class prepended by ``CD_``. 
For example, if you are contributing a class ``MyClass`` the header files for this class must be named :file:`CD_MyClass.H` and the implementation file must be named :file:`CD_MyClass.cpp`.
If your code contains templates or inlined functions, these should be defined in files appended by ``Implem``, e.g. :file:`CD_MyClassImplem.H`.

File headers
------------

Each file shall begin with the following note:

.. code-block:: c++

   /* chombo-discharge
   * Copyright © <Copyright holder 1>
   * Copyright © <Copyright holder 2>     
   * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
   */

where *<Copyright holder 1>*, *<Copyright holder 2>*, etc. are replaced by the copyright holder.

This file header shall be followed by a brief Doxygen documentation, containing at least ``@file``, ``@brief``, and ``@author``.
In addition, include header guards identical to the filename, replacing dots by underscores.
I.e. for a file :file:`CD_MyClass.H` the header guard shall read

.. code-block:: c++

   #ifndef CD_MyClass_H
   #define CD_MyClass_H

   #endif

File inclusions
---------------

File inclusions should use the follow standards for C++, ``Chombo``, and ``chombo-discharge``

1. *C++*. Use brackets, e.g. ``#include <memory>``.
2. ``Chombo``. Use brackets, e.g. ``#include <LevelData.H>``.
3. ``chombo-discharge``. Use brackets and the file name, e.g. ``#include <CD_FieldSolver.H>``.

Example format
--------------

Here is a complete example of a header file in ``chombo-discharge``:

.. code-block:: c++

   /* chombo-discharge
   * Copyright © <Copyright holder 1>
   * Copyright © <Copyright holder 2>     
   * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
   */

   /*!
     @file   CD_MyClass.H
     @brief  This file contains ...
     @author Author name
   */
   
   #ifndef CD_MyClass_H
   #define CD_MyClass_H

   // Std includes (e.g.)
   #include <memory>

   // Chombo includes (e.g.)
   #include <LevelData.H>

   // Our includes (e.g.)
   #include <CD_EBAMRData.H>
   #include <CD_NamespaceHeader.H>

   /*!
     @brief This class does the following: ....
   */
   class MyClass
   {
   public:

   //...
   };

   #include <CD_NamespaceFooter.H>

   #include <CD_MyClassImplem.H> // Inline and template code included at the end. 
   
   #endif

Code syntax
-----------


We use the following syntax:

1. Class names, structs, and namespaces should be in Pascal case where the first letter of every word is capitalized.
   E.g. a class is called ``MyClass``.

2. Class functions should be in Camel case where the first letter of every word but the first is capitalized. 
   E.g. functions should be named ``MyClass::myFunction``

3. Variables should use Pascal-case, with the following requirements:
   
   * Arguments to functions should be prepended by ``a_``. For example ``MyClass::myFunction(int a_inputVariable)``.
     
   * Class members should always be prepended by ``m_``, indicating it is a member of a class. For example ``MyClass::m_functionMember``.
     
   * Static variables are prepended by ``s_``. For example ``MyClass::s_staticFunctionMember``.
     
   * Global variables are prepended by ``//``.

Options files
-------------

Options files are named using the same convention as class files, e.g. ``CD_MyClass.options``.
It is the responsibility of ``MyClass`` to parse these variables correctly.

Everything in the options file should be lower-case, with the exception of the class name which should follow the class name syntax.
If you need a separator for the variable, use an underscore ``_``.
For variables that should be grouped under a common block, one may use a dot ``.`` for grouping them. 
For a class ``MyClass`` and options file might look something like

.. code-block:: text

   MyClass.input_variable = 1.0
   MyClass.bc.x.lo        = dirichlet 1.0

clang-format
------------

We use ``clang-format`` for formatting the source code.
Before opening a pull request for review, navigate to :file:`$DISCHARGE_HOME` and format the code using

.. code-block:: bash

   find Source Physics Geometries Exec \( -name "*.H" -o -name "*.cpp" \) -exec clang-format -i {} +   
