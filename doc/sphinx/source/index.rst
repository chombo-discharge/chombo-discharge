.. PlasmaC documentation master file, created by
   sphinx-quickstart on Thu Jul 26 21:03:53 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

**This is an alpha release of PlasmaC. Development is still in progress, and various bugs may or may not be present. User interfaces can and will change.**

Welcome to `PlasmaC`'s user documentation!
============================================
	      
`PlasmaC` is a scalable computer code for Cartesian two- and three-dimensional simulations of low-temperature fluid advection-diffusion-reaction plasmas in complex geometries. `PlasmaC` features

* Multiphase Poisson and Helmholtz solvers
* Time-dependent multiphase heat solvers
* Diffusive radiative transport
* Convection-diffusion-reaction solvers

The solvers can be run on their own, or they can be coupled through one of our (or your own) physics interfaces. For example, `PlasmaC` can be run as a pure electrostatic solver, a heat or diffusion equation solver, and it is applicable to corona and streamer simulations. 

For scalability, `PlasmaC` is built on top of `Chombo <https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations>`_, and therefore additionally features

* Cut-cell representation of multi-material geometries
* Patch based adaptive mesh refinement
* Good weak and strong scalability to thousands of computer cores

It is our goal that users will be able to use `PlasmaC` without modifying the underlying solvers. Because of this, there are quite general interfaces for describing the plasma physics, setting up boundary conditions, ensuring mesh refinement, and so on. As `PlasmaC` evolves, so will these interfaces. However, we always aim for backward compatibility such that existing `PlasmaC` models can be run on future versions of `PlasmaC`. In particular, there are plans to include Euler equations for the fluid transport and port compute kernels to GPUs. Both of these features will very strongly affect how `PlasmaC` is set up and run. 

This documentation is the user documentation `PlasmaC`. There is also a separate :doxy:`Doxygen API <index>` that can be compiled together with the source code.


.. toctree::
   :maxdepth: 4
   
   Introduction
   Equations
   TemporalDiscretization
   NewSimulations
   Design
   Control
   Tutorials
   ImportantClasses
