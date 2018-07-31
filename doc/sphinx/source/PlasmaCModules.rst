.. _Chap:PlasmaCModules:

PlasmaC modules
---------------

PlasmaC is a modular framework that allows reuse or reimplementation of many classes without affecting the underlying code. In PlasmaC, an entire simulation is run by an object called :ref:`Chap:plasma_engine` . The responsibility of this object is essentially to set up a simulation, run a simulation, and write output files. In fact, the mini-application system in PlasmaC requires the user to construct the plasma_engine object by supplying sufficient supplementary classes through its constructor.

To facilitate code reuse, PlasmaC is set up such that a class :ref:`Chap:plasma_engine` takes as its input a number of modules. In this way, the user can re-use code (for example geometries or physics) over a range of applications. 

.. figure:: figures/plasma_engine.png
   :width: 480px
   :align: center

   PlasmaC modules.

In the above figure, the various modules that go into :ref:`Chap:plasma_engine` are:

* :ref:`Chap:plasma_kinetics` An abstract class that defines the physics. This includes specifying the number of convected species, the number of RTE solvers, and well as specifiying the coupling between all of these. 
* :ref:`Chap:computational_geometry` An implementation of the geometry that will be simulated. Various implementation of this class exist, which you may immediately use. Descriptions of new geometries must be done by the user by either implementing a new computational_geometry class. 
* :ref:`Chap:physical_domain` The physical domain to be simulated. This is a very lightweight class that only describes the axis-aligned box that you wish to simulate. 
* :ref:`Chap:time_stepper` The temporal integrator; the class includes an implementation of the time stepping scheme. The base class is abstract. Implementation of new integration schemes is time consuming, and additional implementations of this class is a developer task. Currently, we support some implicit-explicit schemes and Runge-Kutta schemes. Most users will find the second order Runge-Kutta scheme to be sufficient. This class is the one that owns all the individual solvers. 
* :ref:`Chap:amr_mesh` The AMR mesh engine. This class includes grid generators, coarsening operators, ghost cell interpolation operators and so on. This class is one of the most important base class, and is responsible for orchestratic all spatial operations.
* :ref:`Chap:cell_tagger` (Optional) Class that is responsible for refinement and coarsening decision. The base class is abstract, but users may implement their own classes if they like. 
* :ref:`Chap:geo_coarsener` (Optional) A geometric coarsening class.

You will find a much more thorough explanation of these classes in the :ref:`Chap:ImportantClasses` chapter.
