.. _Chap:Realm:

Realm
=====

``Realm`` is a class for centralizing EBAMR-related grids and operators for a specific AMR grid. 
For example, a ``Realm`` consists of a set of grids (i.e. a ``Vector<DisjointBoxLayout>``) as well as *operators*, i.e., functionality for filling ghost cells or averaging down a solution from a fine level to a coarse level.
One may think of a ``Realm`` as a fully-fledged AMR hierarchy with associated multilevel operators, i.e., how one would usually do AMR.

Dual grid
---------

The reason why ``Realm`` exists at all is due to individual load balancing of algorithmic components. 
The terminology *dual grid* is used when more than one ``Realm`` is used in a simulation, and in this case the user/developer has chosen to solve the equations of motion over a different set of ``DisjointBoxLayout`` on each level.
This approach is very useful when using computational particles since users can load balance the grids for the fluid and particle algorithms separately.
Note that every ``Realm`` consists of the same boxes, i.e., the physical domain and computational grids are the same for all realms. 
The only difference lies primarily in the assignment of MPI ranks to grids, i.e., the load-balancing.

.. _Fig:DualMesh:
.. figure:: /_static/figures/DualMesh.png
   :width: 35%
   :align: center

   Sketch of dual grid approach.
   Each rectangle represents a grid patch and the numbers show MPI ranks. a) Load balancing with the number of grid cells. b) Load balancing with the number of particles.

:numref:`Fig:DualMesh` shows an example of a dual-grid approach.
In  this figure we have a set of grid patches on a particular grid level.
In the top panel the grid patches are load-balanced using the grid patch volume as a proxy for the computational load.
The numbers in each grid patch indicate the MPI rank ownership of the patches.
In the bottom panel we have introduced computational particles in some of the patches.
For particles, the computational load is better defined by the number of computational particles assigned to the patch, and so using the number of particles as a proxy for the load yields different rank ownership over the grid patches.

.. _Chap:RealmHashGrid:

Patch lookup (tile hash grid)
-----------------------------

Each ``Realm`` maintains, on every grid level, a *tile hash grid* that maps a physical position to the grid box (and the owning MPI rank) that covers it.
The hierarchy is decomposed into uniform tiles of ``AmrMesh.min_block_size`` cells, and every grid box is registered under each tile it spans.
Because a box is always a union of aligned ``min_block_size`` tiles -- even when the grids contain variable-sized or anisotropic boxes (see :ref:`Chap:MeshGeneration`) -- the tile-to-box map stays uniform and is stored as a hash map (``std::unordered_map`` keyed by the tile index).
A position is therefore mapped to its owning box in :math:`\mathcal{O}(1)`: the position is converted to a tile by integer division, and that tile is looked up in the per-level hash, finest level first.
A box that spans more than one tile (i.e. larger than ``min_block_size``) is simply registered under each of its tiles; when ``min_block_size`` equals ``max_block_size`` every box is a single tile and the map holds one entry per box.
The hash grid is rebuilt during every regrid, so it always reflects the realm's current grids and rank ownership.

The lookup is exposed through ``Realm::getLevelAndBox``, which returns the finest level whose box covers a given position:

.. code-block:: c++

   const Realm::LevelAndBox result = realm.getLevelAndBox(position);
   // result.level    -- AMR level containing the point (-1 if none)
   // result.boxIndex -- global grid-box index on that level
   // result.rank     -- MPI rank owning that box
   // result.valid    -- false if no box covers the point (i.e. off-domain)

This is the same mapping the particle infrastructure uses to bin particles during ``remap()`` (see :ref:`Chap:Particles`).
``ParticleContainer`` aliases the realm's tile hash grid rather than building its own, so the ``Realm`` is the single source of truth for point-to-patch ownership.

Interacting with realms
-----------------------

Users will not interact with ``Realm`` directly.
Every ``Realm`` is owned by :ref:`Chap:AmrMesh`, and the user will only interact with realms through the public :ref:`Chap:AmrMesh` interface, for example by fetching operators for performing AMR operations. 
It is important, however, to keep track of what data is allocated where.
Fortunately, :ref:`Chap:AmrMesh` will issue plenty of warnings if the user calls a function where the input arguments are realm-wise inconsistent.
