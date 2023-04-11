.. _Chap:Testing:

Code testing
============

To ensure the integrity of ``chombo-discharge``, we include tests in

* :file:`$DISCHARGE_HOME/Exec/Tests` for a functional test suite.
* :file:`$DISCHARGE_HOME/Exec/Convergence` for a code verification tests. 

.. _Chap:TestSuite:

Test suite
----------

To make sure ``chombo-discharge`` compiles and runs as expected, we include a test suite that runs functional tests of many ``chombo-discharge`` components.
The tests are defined in :file:`$DISCHARGE_HOME/Exec/Tests`, and are organized by application type.

Running the test suite
______________________

To do a clean compile and run of all tests, navigate to ``$DISCHARGE_HOME/Exec/Tests`` and execute the following:

.. code-block:: bash

   python3 tests.py --compile --clean --silent --no_exec -cores X

where ``X`` is the number of cores to use when compiling.
By default, this will compile without MPI and HDF5 support.


Advanced configuration
______________________

The following options are available for running the various tests:

* ``--compile`` Compile all tests. 
* ``--clean`` Do a clean recompilation.
* ``--silent`` Turn off terminal output.
* ``--benchmark`` Generate benchmark files.
* ``--no_exec`` Compile, but do not run the test.
* ``--compare`` Run, and compare with benchmark files.
* ``-dim <number>`` Run only the specified 2D or 3D tests.  
* ``-mpi <true/false>`` Use MPI or not.
* ``-hdf <true/false>`` Use HDF5 or not.  
* ``-cores <number>`` Run with specified number of cores.
* ``-suites <string>`` Run a specific application test suite.
* ``-tests <string`` Run a specific test.

For example, to compile with MPI and HDF5 enabled:

.. code-block:: text

   python3 tests.py --compile --clean --silent --no_exec -mpi=true -hdf=true -cores 8

If one only wants to compile a specified test suite:

.. code-block:: text

   python3 tests.py --compile --clean --silent --no_exec -mpi=true -hdf=true -cores 8 -suites Electrostatics

One can also restrict the tests to specific dimensionality, e.g.

.. code-block:: text

   python3 tests.py --compile --clean --silent --no_exec -mpi=true -hdf=true -cores 8 -suites Electrostatics -dim=2

Using benchmark files
_____________________

The test suite can generate benchmark files which can later be compared against new test suite output files.
This is often a good idea if one wants to ensure that changes the ``chombo-discharge`` code does not unintentionally change simulations. 
In this case one can run the test suite and generate benchmark files *before* adding changes to ``chombo-discharge``.
Once the code development is completed, the benchmark files can later be bit-wise (using `h5diff <https://support.hdfgroup.org/HDF5/doc/RM/Tools/h5diff.htm>`_) compared against the results of a later test suite.

This consists of the following steps:

#. *Before* making changes to ``chombo-discharge``, generate benchmark files with

   .. code-block:: bash

      python3 tests.py --compile --clean --silent --benchmark -mpi=true -hdf=true -cores X		   

#. Make the required changes to the ``chombo-discharge`` code.

#. Run the test suite again, and compare benchmark and output files as follows:

   .. code-block:: bash

      python3 tests.py --compile --clean --silent --compare -mpi=true -hdf=true -cores X		   

When running the tests this way, the output files are bit-wise compared and a warning is issued if the files not exactly match. 

.. _Chap:AutomatedTests:      

Automated testing
-----------------

On `GitHub <https://github.com/chombo-discharge/chombo-discharge>`_, the test suite is integrated with GitHub actions and are automatically run when opening a pull request for review. 
In general, all tests must pass before a pull request can be merged.
The test status can be observed either in the pull request, or at `<https://github.com/chombo-discharge/chombo-discharge/actions>`_.
The automated tests run ``chombo-discharge`` with ``DEBUG=TRUE`` and ``OPT=FALSE`` in order to catch assertion errors or other code breaks.
They usually take 1-2 hours to complete.

The automated tests will clone, build, and run the ``chombo-discharge`` test suite for various configurations:

* Parallel and serial.
* With or without HDF5.
* In 2D and 3D.

The tests are run with the following compiler suites:

* GNU.
* Intel oneAPI.

.. _Chap:ConvergenceTests:  

Convergence testing
-------------------

To ensure that the various components in ``chombo-discharge`` converge at desired truncation order, many modules are equipped with their own convergence tests.
These are located in :file:`$DISCHARGE_HOME/Exec/Convergence`.
The tests are too extensive to include in continuous integration, and they must be run locally like a regular ``chombo-discharge`` application.
Our approach for convergence testing is found in :ref:`Chap:VV`.
