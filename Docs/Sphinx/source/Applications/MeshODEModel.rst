.. _Chap:MeshODEModel:

Mesh ODE
========

The ``MeshODEStepper`` module is designed to illustrate the usage of a :ref:`Chap:MeshODESolver`.
This model has no functionality for cell-refinement based on the solutions, so refinement is limited to the geometric refinement methods supplied by :ref:`Chap:Driver`.

The main class implements :ref:`Chap:TimeStepper` and is templated as follows:

.. literalinclude:: ../../../../Physics/MeshODE/CD_MeshODEStepper.H
   :language: c++
   :lines: 36-40
   :dedent: 4

Here, ``N`` is the number of variables that will be stored on the mesh.

Example problem
---------------

The example problem set up by ``MeshODEStepper`` is

.. math::

   \partial_t\phi(t) = \cos\left(2\pi f t\right),

where :math:`f` is a user-supplied frequency.
The user can also set the initial value for :math:`\phi` through the configuration options for ``MeshODEStepper<N>``, see :ref:`Chap:MeshODEStepperConfiguration`

Setting up a new problem
------------------------

To set up a new problem, using the Python setup tools in :file:`$DISCHARGE_HOME/Physics/MeshODEStepper` is the simplest way.
A full description is available in the ``README.md`` contained in the folder:

.. literalinclude:: ../../../../Physics/MeshODEStepper/README.md
   :language: markdown
	      
To see available setup options, use

.. code-block:: bash

   ./setup.py --help

.. _Chap:MeshODEStepperConfiguration:

Solver configuration
--------------------

The ``MeshODEStepper`` class comes with configurable input options that can be adjusted at runtime, which are listed below

.. literalinclude:: ../../../../Physics/MeshODE/CD_MeshODEStepper.options
   :language: text   
