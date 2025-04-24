.. _Chap:AmrMesh:

AmrMesh
========

:ref:`Chap:AmrMesh` handles (almost) all spatial operations in ``chombo-discharge``.
Internally, :ref:`Chap:AmrMesh` contains a bunch of operators that are useful across classes, such as ghost cell interpolation operators, coarsening operators, and stencils for interpolation and extrapolation near the embedded boundaries.
:ref:`Chap:AmrMesh` also contains routines for generation and load-balancing of grids based and also contains simple routines for allocation and deallocation of memory.

.. note::
   :ref:`Chap:AmrMesh` only handles spatial *operations*, it otherwise has limited knowledge of numerical discretizations. 

:ref:`Chap:AmrMesh` is an integral part of ``chombo-discharge``, and users will never have the need to modify it unless they are implementing something entirely new.
The behavior of :ref:`Chap:AmrMesh` is modified through its available input parameters.

.. tip::

   `AmrMesh C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classAmrMesh.html>`_

Main functionality
------------------

There are two main functionalities in ``AmrMesh``:

1. Building grid hierarchies, and providing geometric information
2. Providing AMR operators.

The grids in ``AmrMesh`` consist of a ``DisjointBoxLayout`` (see :ref:`Chap:Basics`) on each level, supported also by the EB information (``EBISLayout``).
Recall that each grid patch in a ``DisjointBoxLayout`` is owned by a unique rank.
However, since ``chombo-discharge`` supports multiple decompositions, we support the use of multiple ``DisjointBoxLayout`` describing the same grid level.
Although these ``DisjointBoxLayout`` consist of the same boxes, the patch-to-rank mapping can be different (see :ref:`Chap:Realm`).
To fetch a grid on a particular level, one can call ``AmrMesh::getGrids(const std::string a_realm)``.
E.g.

.. code-block:: c++

   const std::string myRealm;		
   const int myLevel;
   
   const DisjointBoxLayout& dbl = m_amr->getGrids("myRealm")[myLevel];

Likewise, to fetch the geometric (EB) information for a specified realm, phase, and level:

.. code-block:: c++

   const std::string myRealm;
   const int myLevel;
   const phase::which_phase myPhase;
   
   const EBISLayout& ebisl = m_amr->getEBISLayout(myRealm, myPhase)[myLevel];

   
In addition to the grids, the user can fetch AMR operators.
These are, for example, coarsening operators, interpolation operators, ghost cell interpolators etc.
To save some regrid time, we don't always build every AMR operator that we might ever need, but have solvers *register* the ones that they specifically need.
See :ref:`Chap:Realm` for details.

The API of ``AmrMesh`` is quite extensive, and it has functionality both for providing grid operators as well as:

* Maintaining overview of grids and operators.
* Allocating grid and particle data.
* Interpolating to new grids.
* Coarsening data.
* Updating ghost cells.
* Interpolating data to e.g. EB centroids.

The `AmrMesh API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classAmrMesh.html>`_ is a good place to start for figuring out the ``AmrMesh`` functionality.

.. _Chap:AmrParticleIntersection:

Particle intersection
---------------------

Class options
-------------

The class options below control ``AmrMesh``:

* ``AmrMesh.lo_corner``. Low corner of problem domain (e.g. 0 0 0)
* ``AmrMesh.hi_corner``. High corner of problem domain (e.g. 1 1 1). 
* ``AmrMesh.verbosity``. Class verbosity. Leave to -1 unless you are debugging. 
* ``AmrMesh.coarsest_domain``. Number of grid cells on coarsest domain
* ``AmrMesh.max_amr_depth``. Maximum number of refinement levels. 
* ``AmrMesh.max_sim_depth``. Maximum simulation depth.
  Values :math:`< 0` means that grids can be generated with depths up to ``AmrMesh.max_amr_depth``. 
* ``AmrMesh.fill_ratio``. Fill ratio for BR grid generation
* ``AmrMesh.buffer_size``. Buffer size for BR grid generation. 
* ``AmrMesh.grid_algorithm``. Grid generation algorithm. Valid options are *br* or *tiled*. See :ref:`Chap:MeshGeneration` for details. 
* ``AmrMesh.box_sorting``. Box sorting algorithm. Valid options are *std*, *morton*, or *shuffle*. 
* ``AmrMesh.blocking_factor``. Blocking factor. 
* ``AmrMesh.max_box_size``. Maximum box size. 
* ``AmrMesh.max_ebis_box``. Maximum box size during EB geometry generation. 
* ``AmrMesh.ref_rat``. Refinement ratios. 
* ``AmrMesh.num_ghost``. Number of ghost cells for mesh data. 
* ``AmrMesh.lsf_ghost``. Number of ghost cells when allocating level-set function on the grid. 
* ``AmrMesh.eb_ghost``. Number of ghost cells for EB moments. 
* ``AmrMesh.centroid_interp``. Which centroid interpolation stencils to use. Good options are *minmod*, *linear*.
* ``AmrMesh.eb_interp``. EB interpolation stencils. Good options are *minmod*, *pwl*
* ``AmrMesh.redist_radius``. Redistribution radius. 

.. warning::

   ``chombo-discharge`` only supports uniform resolution (i.e., cubic grid cells).
   I.e. the user must specify consist domain sizes and resolutions.   

Runtime options
---------------

The following options are runtime options for ``AmrMesh``: 

* ``AmrMesh.verbosity``. 
* ``AmrMesh.fill_ratio``. 
* ``AmrMesh.irreg_growth``. 
* ``AmrMesh.buffer_size``. 
* ``AmrMesh.grid_algorithm``. 
* ``AmrMesh.box_sorting``. 
* ``AmrMesh.blocking_factor``. 
* ``AmrMesh.max_box_size``.
* ``AmrMesh.centroid_interp``
* ``AmrMesh.eb_interp``  

These options only affect the grid generation method and parameters, and are thus only effective after the next regrid.


