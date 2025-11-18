.. _Chap:GeometryModel:

Geometry
========

Functionality
-------------

The ``GeometryStepper`` class is a solver-less class used for supporting the development of new geometries, which in general derive from :ref:`Chap:ComputationalGeometry`.
``GeometryStepper`` implements :ref:`Chap:TimeStepper` without any true functionality.
E.g., calling the ``advance`` method simply returns an infinitely large time step, and ``GeometryStepper`` provides no I/O functionality or cell refinemention functionality.

.. important::

   Using ``GeometryStepper`` is the simplest way of creating new geometries since all other functionality other than geometry generation is essentially turned off.

Setting up a new problem
------------------------

To set up a new problem, using the Python setup tools in :file:`$DISCHARGE_HOME/Physics/Geometry` is the simplest way.
A full description is available in the ``README.md`` contained in the folder:

.. literalinclude:: ../../../../Physics/Geometry/README.md
   :language: markdown
	      
To see available setup options, use

.. code-block:: bash

   python setup.py --help
