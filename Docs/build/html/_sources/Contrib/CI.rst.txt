.. _Chap:CI:

Continuous integration
======================

``chombo-discharge`` uses continuous integration at GitHub.
These tests consist of cloning, installing, and running various ``chombo-discharge`` tests.

GitHub actions
--------------

The tests defined in ``$DISCHARGE_HOME/Exec/Tests`` are automatically run when opening a pull request.
In general, all tests must pass before a pull request can be merged.
The tests status can be observed either in the pull request, or at `<https://github.com/chombo-discharge/chombo-discharge/actions>`_.
GitHub actions usually take 1-2 hours to complete.

Running tests locally
---------------------

The tests are defined in ``$DISCHARGE_HOME/Exec/Tests`` and can be run using various configurations.
Generally, it is a good idea to run the tests locally before running them at GitHub (which can be slower).

Running all tests
_________________

To do a clean compile and run of all tests, navigate to ``$DISCHARGE_HOME/Exec/Tests`` and execute the following:

.. code-block:: bash

   python3 tests.py --compile --clean --silent --no_compare --parallel -cores X

where ``X`` is the number of cores to use when compiling and running.
If the tests should be run in serial, omit the ``--parallel -cores X`` flag.

Test suite options
__________________

The following options are available for running the various tests:

Using benchmark files
_____________________

The test suite can generate benchmark files which can later be compared against new test suite output files.
This can be a good idea if one wants to ensure that changes the ``chombo-discharge`` code does not affect output files.
In this case one can run the test suite and generate benchmark files *before* adding changes to ``chombo-discharge``.
Once the code development is completed, the benchmark files can later be bit-wise compared against the results of a later test suite.

This consists of the following steps:

#. *Before* making changes to ``chombo-discharge``, generate benchmark files with

   .. code-block:: bash

      python3 tests.py --compile --clean --silent --parallel -cores X --benchmark

#. Make the required changes to the ``chombo-discharge`` code.

#. Run the test suite again, and compare benchmark and output files as follows:

   .. code-block:: bash

      python3 tests.py --compile --clean --silent --parallel -cores X
