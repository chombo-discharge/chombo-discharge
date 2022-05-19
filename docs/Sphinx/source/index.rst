.. chombo-discharge documentation master file, created by
   sphinx-quickstart on Thu Jul 26 21:03:53 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. important::

   **This is an alpha release. Development is still in progress, and various bugs may be present. User interfaces can and will change.**

Welcome to ``chombo-discharge``'s user documentation
****************************************************
	      
``chombo-discharge`` is a modular and scalable research code for Cartesian two- and three-dimensional simulations of low-temperature plasmas in complex geometries.
The code is hosted at `GitHub <https://github.com/chombo-discharge/chombo-discharge>`_ together with the source files for this documentation.

``chombo-discharge`` supports

* Fully written in C++.
* Support for complex geometries.
* Scalar advection-diffusion-reaction processes.
* Electrostatics with support for electrodes and dielectrics.
* Radiative transport as a diffusion or Monte Carlo process.
* Particle-mesh operations (like Particle-In-Cell)
* Parallel I/O with HDF5.
* Dual-grid simulations with individual load balancing of fluid and particles.
* Various multi-physics interfaces that use the above solvers.
* Various time integration schemes.  

.. * Îto particle models for microscopic drift-diffusion-reaction processes.  

Numerical solvers are designed to run either on their own, or as a part of a larger application. 

For scalability, ``chombo-discharge`` is built on top of `Chombo 3 <https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations>`_, and therefore additionally features

* Cut-cell representation of multi-material geometries.
* Patch based adaptive mesh refinement.
* Weak and strong scalability to thousands of computer cores.

Our goal is that users will be able to use ``chombo-discharge`` without modifying the underlying solvers.
There are interfaces for describing e.g. the plasma physics, boundary conditions, mesh refinement, etc.
As ``chombo-discharge`` evolves, so will these interfaces.
We aim for (but cannot guarantee) backward compatibility such that existing ``chombo-discharge`` models can be run on future versions of ``chombo-discharge``.

.. This is for getting rid of the TOC in html view. 
.. raw:: html

   <style>
   /* front page: hide chapter titles
    * needed for consistent HTML-PDF-EPUB chapters
    */
   div#introduction.section,
   div#discretization.section,
   div#design.section,
   div#solvers.section,
   div#physics-models.section,
   div#tutorial.section,
   div#utilities.section,
   div#contributing.section,
   div#references.section,            
   div#epilogue.section {
       display:none;
   }
   </style>

Introduction
*************

.. toctree::
   :maxdepth: 3
   :caption: Introduction
   :hidden:

   Base/Documentation
   Base/Overview
   Base/GettingStarted
   Base/Visualization
   Base/Control

Discretization
**************

.. toctree::
   :maxdepth: 3
   :caption: Discretization
   :hidden:

   Source/SpatialDiscretization      
   Source/ChomboBasics      
   Source/MeshData
   Source/Particles
   Source/Realm
   Source/LinearSolvers

Physics models
**************

.. toctree::
   :maxdepth: 3
   :caption: Applications
   :hidden:	     

   Applications/ImplementedModels
   Applications/AdvectionDiffusionModel
   Applications/BrownianWalkerModel         
   Applications/CdrPlasmaModel
   Applications/ElectrostaticsModel
   Applications/GeometryModel      
   Applications/RadiativeTransferModel      

..   ItoPlasmaModel

Tutorial
********

.. toctree::
   :maxdepth: 3
   :caption: Tutorial
   :hidden:	     

   Tutorials/Tutorial   

Design
******

.. toctree::
   :maxdepth: 3
   :caption: Design
   :hidden:	     

   Source/Driver
   Source/ComputationalGeometry
   Source/TimeStepper
   Source/AmrMesh   
   Source/CellTagger
   Source/GeoCoarsener

Solvers
*******

.. toctree::
   :maxdepth: 3
   :caption: Solvers
   :hidden:	     

   Solvers/CDR
   Solvers/Electrostatics
   Solvers/RTE
   Solvers/Sigma
   Solvers/Ito

Utilities
*********  

.. toctree::
   :maxdepth: 3
   :caption: Utilities
   :hidden:	     

   Utilities/LookupTable
   Utilities/RandomNumbers
   Utilities/LeastSquares
   Utilities/ComplexGeometries

Contributing
************  

.. toctree::
   :maxdepth: 3
   :caption: Contributing
   :hidden:	     

   Contrib/Contributions
   Contrib/CodeStandard

References
**********

.. toctree::
   :maxdepth: 3
   :caption: References
   :hidden:	     

   ZZReferences
