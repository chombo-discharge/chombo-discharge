.. PlasmaC documentation master file, created by
   sphinx-quickstart on Thu Jul 26 21:03:53 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PlasmaC's user documentation!
========================================

PlasmaC is a scalable computer code for two- and three-dimensional simulations of low-temperature fluid plasmas in complex geometries. PlasmaC supports

* Multi-material Poisson solvers
* Diffusion radiation transpor tsolvers
* Convection diffusion reaction systems

PlasmaC is built on top of Chombo, and therefore additionally features

* Cut-cell representation of multi-material geometries
* Patch based adaptive mesh refinement
* Excellent scalability from 1 to 100,000 cores.

Take a look at the gallery for some examples. 

This documentation is the user documentation for usage and examples on how to use PlasmaC. There is a separate doxygen API documentation of the source code. 

.. toctree::
   :maxdepth: 1

   Introduction
   WorkedExample
   Interface
   CodeStructure
   ImportantClasses
   Gallery

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
