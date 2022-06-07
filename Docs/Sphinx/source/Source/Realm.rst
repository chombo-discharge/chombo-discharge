.. _Chap:Realm:

Realm
=====

``Realm`` is a class for centralizing EBAMR-related grids and operators for a specific AMR grid. 
For example, a ``Realm`` consists of a set of grids (i.e. a ``Vector<DisjointBoxLayout>``) as well as *operators*, e.g. functionality for filling ghost cells or averaging down a solution from a fine level to a coarse level.
One may think of a ``Realm`` as a fully-fledged AMR hierarchy with associated multilevel operators, i.e. how one would usually do AMR.



Dual grid
---------

The reason why ``Realm`` exists at all is due to individual load balancing of algorithmic components. 
The terminology *dual grid* is used when more than one ``Realm`` is used in a simulation, and in this case the user/developer has chosen to solve the equations of motion over a different set of ``DisjointBoxLayout`` on each level.
This approach is very useful when using computational particles since users can quickly generate separate Eulerian sets of grids for fluids and particles, and the grids can then be load balanced separately.
Note that every ``Realm`` consists of the same boxes, i.e. the physical domain and computational grids are the same for all realms. 
The difference lies primarily in the assignment of MPI ranks to grids; i.e. the load-balancing and domain decomposition.

.. _Fig:DualMesh:
.. figure:: /_static/figures/DualMesh.png
   :width: 35%
   :align: center

   Sketch of dual grid approach.
   Each rectangle represents a grid patch and the numbers show MPI ranks. a) Load balancing with the number of grid cells. b) Load balancing with the number of particles.

:numref:`Fig:DualMesh` shows an example of a dual-grid approach.
In  this figure we have a set of grid patches on a particular grid level.
In the top panel the grid patches are load-balanced using the grid patch volume as a proxy for the computational load.
The numbers in each grid patch indicates the MPI rank ownership of the patches.
In the bottom panel we have introduced computational particles in some of the patches.
For particles, the computational load is better defined by the number of computational particles assigned to the patch, and so using the number of particles as a proxy for the load yields different rank ownership over the grid patches.

Interacting with realms
-----------------------

Users will not interact with ``Realm`` directly.
Every ``Realm`` is owned by ``AmrMesh``, and the user will only interact with realms through the public ``AmrMesh`` interface, for example by fetching operators for performing AMR operations. 
In addition, data that is defined on one realm can be copied to another; ``EBAMRData<T>`` takes care of this.
You will simply call a copier:

.. code-block:: c++

   EBAMRCellData realmOneData;
   EBAMRCellData realmTwoData;

   realmOneData.copy(realmTwoData);

The rest of the functionality uses the public interface of :ref:`Chap:AmrMesh`.
For example for coarsening of multifluid data:

.. code-block:: c++

   std::string multifluidRealm;
   MFAMRCellData multifluidData;
   AmrMesh amrMesh;

   amrMesh.averageDown(multifluidData, multifluidRealm);
