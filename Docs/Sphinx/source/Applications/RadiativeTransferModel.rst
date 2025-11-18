.. _Chap:RadiativeTransferModel:

Radiative transfer
==================

The radiative transfer model runs a single solver using the :ref:`Chap:RtSolver` interface, where the underlying solver can derived from any subclass of :ref:`Chap:RtSolver`.
This module is designed to show how to set up and run :ref:`Chap:RtSolver`, and the code was not written to solve any particular problem.
Selecting between different types of solvers must be done at compile time. 

.. tip::

   The source code is located in :file:`$DISCHARGE_HOME/Physics/RadiativeTransfer`.

The model consists of the following implementation files:

* :file:`CD_RadiativeTransferStepper.H` which implements :ref:`Chap:TimeStepper`. 
* :file:`CD_RadiativeTransferSpecies.H` which implements the necessary transport conditions through :ref:`Chap:RtSpecies`. 

.. note::

   The current radiative transfer module does not incorporate solver-based adaptive mesh refinement, so refinement is restricted to refinement and coarsening through the :ref:`Chap:Driver` interface.

Basic problem
-------------

Currently, ``RadiativeTransferStepper`` simply instantiates a solver and advances the solution using a synthetic source given by

.. math::

   \eta\left(\mathbf{x}\right) = \eta_0\exp\left[-\frac{\left(\mathbf{x}-\mathbf{x}_0\right)^2}{2 R^2}\right],

where the source strength :math:`\eta_0`, source radius :math:`R`, and source center :math:`\mathbf{x}_0` are set by the user.
Similarly, the user can set the photon absorption length, see :ref:`Chap:RadiativeTransferStepperConfiguration`.

.. _Chap:RadiativeTransferStepperConfiguration:

Solver configuration
--------------------

The ``RadiativeTransferStepper`` class comes with user-configurable input options that can be adjusted at runtime.
These configuration options are given below.

.. literalinclude:: ../../../../Physics/RadiativeTransfer/CD_RadiativeTransferStepper.options
   :language: text

Setting up a new problem
------------------------

To set up a new problem, using the Python setup tools in :file:`$DISCHARGE_HOME/Physics/RadiativeTransfer` is the simplest way.
A full description is available in the ``README.md`` contained in the folder:

.. literalinclude:: ../../../../Physics/RadiativeTransfer/README.md
   :language: markdown
	      
To see available setup options, use

.. code-block:: bash

   python setup.py --help


.. tip::

   To set up a new problem using either a Monte-Carlo or diffusion-based radiative transfer solver, pass in the flag ``-RtSolver`` to set up the problem using the desired solver. 
