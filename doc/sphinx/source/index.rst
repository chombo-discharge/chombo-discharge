.. PlasmaC documentation master file, created by
   sphinx-quickstart on Thu Jul 26 21:03:53 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to PlasmaC's user documentation!
========================================
	      
PlasmaC is a scalable computer code for two- and three-dimensional simulations of low-temperature fluid plasmas in complex geometries. PlasmaC features

* Multi-material Poisson solvers
* Diffusion radiation transport
* Convection diffusion reaction systems

PlasmaC is built on top of `Chombo <https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations>`_, and therefore additionally features

* Cut-cell representation of multi-material geometries
* Patch based adaptive mesh refinement
* Excellent scalability from 1 to 100,000 cores.

Take a look at the :ref:`Chap:Gallery` for some examples. 

This documentation is the user documentation for usage and examples on how to use PlasmaC. There is also a separate :doxy:`Doxygen API <index>` documentation of the source code. 

Table of contents
=================

.. toctree::
   :maxdepth: 3

   Obtaining
   Introduction
   Interface
   CreatingGeometries
   NewSimulations
   RunningSimulations
   ImportantClasses
   Solvers
   FAQ
   Bugs
   Todo
   Visualization
   Gallery

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _Chap:contact:

Contact us
==========

Send us an email at robert.marskar@sintef.no
