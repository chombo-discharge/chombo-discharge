.. _Chap:Contributions:

Contributions
=============

We welcome feedback, bug reports, and contributions to ``chombo-discharge``.
If you have feedback, questions, or general types of queries, use the issue tracker or discussion tab at `https://github.com/chombo-discharge <https://github.com/chombo-discharge/chombo-discharge>`_.

Pull requests
-------------

If you want to submit code to ``chombo-discharge``, use the pull request system at `https://github.com/chombo-discharge <https://github.com/chombo-discharge/chombo-discharge>`_.
When submitting a pull request, mark it as *draft* if you believe the pull request is not yet ready for merging.
Note that when a pull request is submitted, you are asked to provide a brief summary of the changes that are made.

.. important::

   Always squash your git commits into a single commit for the PR as a whole.

It will be beneficial to first compile the changes locally and run e.g. the test suite, see :ref:`Chap:Testing` in debug mode.
Also make sure to run e.g. *valgrind* to spot memory leaks.
In serial, running valgrind can be done by calling valgrind with ``--leak-check=yes``, e.g.:

.. code-block:: bash

   valgrind --leak-check=full --track-origins=yes <my_executable> <my_inputs_file>		

With MPI this looks like

.. code-block:: bash

   mpirun -np 12 valgrind --leak-check=full --track-origins=yes <my_executable> <my_inputs_file>		

Bug reports
-----------

``chombo-discharge`` is probably not bug-free.
If encountering unexpected behavior, do not hesitate to use the issue tracker at `<https://github.com/chombo-discharge/chombo-discharge/issues>`_.

Continuous integration
----------------------

``chombo-discharge`` uses continuous integration (CI) with GitHub actions for:

* Running the test suite (see :ref:`Chap:Testing`). 
* Building HTML, PDF, and doxygen documentation.
* Ensuring correct code format.

When submitting pull request for review, the above tests will start.
In general, all the above tests should pass before merging the pull request into the main branch.
Note that all actions on GitHub run in debug mode (``DEBUG=TRUE``) for turning on assertions.


Upon merging with the main branch, the documentation is again rebuilt and deployed to GitHub pages (see `GitHub pages <https://pages.github.com/>`_).
This ensures that the online documentation (HTML, PDF, and doxygen) is always up-to-date with the latest ``chombo-discharge`` release.

