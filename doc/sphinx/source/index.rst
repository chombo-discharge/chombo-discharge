.. PlasmaC documentation master file, created by
   sphinx-quickstart on Thu Jul 26 21:03:53 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

**This is an alpha release of PlasmaC. Development is still in progress, and various bugs may or may not be present. User interfaces can and will change.**

Welcome to `PlasmaC`'s user documentation!
============================================
	      
`PlasmaC` is a modular and scalable computer code for Cartesian two- and three-dimensional simulations of low-temperature plasmas in complex geometries. `PlasmaC` evolves deterministic or fluctuating equations of hydrodynamics and features

* Electrostatic field solver
* Stationary or transient diffusive RTE
* Stationary or transient Monte Carlo based RTE
* Deterministic and fluctuating advection-diffusion-reaction solvers
* Various time integration schemes
    
   * Godunov splitting
   * Strang splitting based on high-order SSPRK schemes
   * Semi-implicit spectral deferred corrections (SISDC)
   * Modified Euler-Maruyama for fluctuating hydrodynamics
* Parallel I/O with HDF5
* Sensible and simple-to-use physics interfaces

Solvers can be run on their own, or they can be coupled through our physics interfaces. 

For scalability, `PlasmaC` is built on top of `Chombo <https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations>`_, and therefore additionally features

* Cut-cell representation of multi-material geometries
* Patch based adaptive mesh refinement
* Good weak and strong scalability to thousands of computer cores

Our goal is that users will be able to use `PlasmaC` without modifying the underlying solvers. There are interfaces for describing the plasma physics, setting up boundary conditions, ensuring mesh refinement, and so on. As `PlasmaC` evolves, so will these interfaces. We aim for (but cannot guarantee) backward compatibility such that existing `PlasmaC` models can be run on future versions. 

This documentation is the user documentation `PlasmaC`. There is also a separate :doxy:`Doxygen API <index>` that can be compiled together with the source code.


.. toctree::
   :maxdepth: 4
   
   GettingStarted
   Model
   NewSimulations
   Design
   Control
   TemporalDiscretization
   Tutorials
   ImportantClasses
