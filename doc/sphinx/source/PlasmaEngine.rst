.. _Chap:plasma_engine:

plasma_engine
-------------

:ref:`Chap:plasma_engine` is the most important class in PlasmaC, and is the object that choreographs a mini-app. In order to be instantiated, :ref:`Chap:plasma_kinetics` must be fed all the other modules that we discuss below. On the whole, :ref:`Chap:plasma_engine` is given the following responsibilities:

* Given a geometry, it will call for generation of the necessary geometric information for the EB grid. This is done by fetching the geometry from :ref:`Chap:computational_geometry`, the physical simulation domain from :ref:`Chap:physical_domain`, and the grid information from :ref:`Chap:amr_mesh`. 
* It is reponsible for setting up simulations for fresh starts, or for restarts (by reading a checkpoint file). See :ref:`Chap:RestartingSimulations` for details. 
* It calls for instantiation of solvers at appropriate times. This is done by calls to :ref:`Chap:time_stepper`, which has the ownership of the solvers.
* It keeps track of the time step, and ensures that the time in all solvers are synchronized.
* The class is responsible for regridding, which occurs through two separate mechanisms: Tagging of irregular cells and tagging of regular cells. The former is only done once, while the latter is controlled through an input parameter. For regular cells, :ref:`Chap:plasma_engine` will perform calls to (derived classes of) :ref:`Chap:cell_tagger`.
* :ref:`Chap:plasma_engine` does all the output and checkpointing.

Here are the options for plasma_engine in the current version of the code:

.. literalinclude:: links/plasma_engine.options

Most options here are self-explanatory. However, we explicitly mention a few that may not be immediately clear. Firstly, **geometry_only** is a special option that *only* generates the geometry. It will be written to an HDF5 file whose name depends on **output_directory** and **output_names**. When you develop your applications, it is often convenient to set this one to *true*, since this will skip a bunch of initialization stages (such as solving the Poisson equation, for example). It allows for quick debugging of your geometry.

The **restart** flag allows you to restart a simulation from a certain checkpoint step. See :ref:`Chap:RestartingSimulations` for details.

The geometry refinement options specify to which level the geometry will be refined. Internally, this is done by tagging irregular cells down to certain levels. **refine_geometry** is a master option, refining all surfaces down to the specified level. The remaining options allow you to individually tag certain surfaces, such as dielectric-gas surfaces, electrode-gas surfaces, electrode-dielectric surfaces and so on. We remark that :ref:`Chap:geo_coarsener` offers a way of leveraging this coarse-grained refinement by removing tags in certain places. 
