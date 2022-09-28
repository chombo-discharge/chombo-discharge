.. _Chap:Particles:

Particles
=========

``chombo-discharge`` supports computational particles using native ``Chombo`` particle data.
The source code for the particle functionality resides in :file:`$DISCHARGE_HOME/Source/Particle`.

.. _Chap:GenericParticle:

GenericParticle
---------------

.. _Chap:ParticleContainer:

ParticleContainer
------------------

The ``ParticleContainer<P>`` is a template class that

#. Stores computational particles of type ``P`` over an AMR hierchy.
#. Provides infrastructure for mapping and remapping. 

``ParticleContainer<P>`` uses the ``Chombo`` structure ``ParticleData<P>`` under the hood, and therefore has template constraints on ``P``.
The simplest way to use ``ParticleContainer`` for a new type of particle is to let ``P`` inherit from the ``Chombo`` class ``BinItem``.
``BinItem`` only has a single member variable which is its position, but derived classes will contain more and must therefore also add new linearization functions if the new member variables should be communicated.
There are many examples of ``chombo-discharge`` particles, see e.g. `TracerParticle <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classTracerParticle.html>`_ or `Photon <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classPhoton.html>`_.
Please refer to the ``Chombo`` design document for complete specification on the template constraints of ``P``, or see some of the examples in ``chombo-discharge``. 

Data structures
---------------

List<P> and ListBox<P>
______________________

At the lowest level the particles are always stored in a linked list ``List<P>``.
The class can be simply be through of as a regular list of ``P`` with non-random access. 

The ``ListBox<P>`` consists of a ``List<P>`` *and* a ``Box``.
The latter specifies the grid patch that the particles are assigned to.

To get the list of particles from a ``ListBox<P>``:

.. code-block::

   ListBox<P> myListBox;
   
   List<P>& myList = myListBox.listItems();


ListIterator<P>
_______________

In order to iterate over particles, use an iterator ``ListIterator<P>`` (which is not random access):

.. code-block:: c++

   List<P> myParticles;
   for (ListIterator<P> lit(myParticles); lit.ok(); ++lit){
      P& p = lit();
      
      // ... do something with this particle
   }

ParticleData<P>
_______________

On each grid level, ``ParticleContainer<P>`` stores the particles in a ``Chombo`` class ``ParticleData``. 

.. code-block:: c++

   template <class P>
   ParticleData<P>

where ``P`` is the particle type.
``ParticleData<P>`` can be thought of as a ``LevelData<ListBox<P> >``, although it actually inherits from ``LayoutData<ListBox<P> >``.
Each grid patch contains a ``ListBox<P>`` of particles. 


AMRParticles<P>
_______________

``AMRParticles<P>`` is our AMR version of ``ParticleData<P>``.
It is a simply a typedef of a vector of pointers to ``ParticleData<P>`` on each level:

.. code-block:: c++

   template <class P>
   using AMRParticles = Vector<RefCountedPtr<ParticleData<P> > >;

Again, the ``Vector`` indicates the AMR level and the ``ParticleData<P>`` is a distributed data holder that holds the particles on each AMR level.

Basic use
---------

Here, we give some examples of basic use of ``ParticleContainer``.
For the full API, see the `ParticleContainer doxygen documentation <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classParticleContainer.html>`_.

Getting the particles
_____________________

To get the particles from a ``ParticleContainer<P>`` one can call ``AMRParticles<P>& ParticleContainer<P>::getParticles()`` which will provide the particles:

.. code-block:: c++

   ParticleContainer<P> myParticleContainer;
   
   AMRParticles<P>& myParticles = myParticleContainer.getParticles();

Alternatively, one can fetch directly from a specified grid level as follows:

.. code-block:: c++

   int lvl;
   ParticleContainer<P> myParticleContainer;
   
   ParticleData<P>& levelParticles = myParticleContainer[lvl];

Iterating over particles
________________________

To do something basic with the particle in a ``ParticleContainer<P>``, one will typically iterate over the particles in all grid levels and patches.

The code bit below shows a typical example of how the particles can be moved, and then remapped onto the correct grid patches and ranks if they fall off their original one. 

.. code-block:: c++

   ParticleContainer<P> myParticleContainer;

   // Iterate over grid levels
   for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){

      // Get the grid on this level. 
      const DisjointBoxLayout& dbl = m_amr->getGrids(myParticleContainer.getRealm())[lvl];

      // Get the distributed particles on this level
      ParticleData<P>& levelParticles = myParticleContainer[lvl]

      // Iterate over grid patches on this level
      for (DataIterator dit(dbl); dit.ok(); ++dit){

         // Get the particles in the current patch.
	 List<P>& patchParticles = levelParticles[dit()].listItems();

	 // Iterate over the particles in the current patch.
	 for (ListIterator<P> lit(patchParticles); lit.ok(); ++lit){
	    P& p = lit();

	    // Move the particle
	    p.position() = ...
	 }
      }
   }

   // Remap particles onto new patches and ranks (they may have moved off their original ones)
   myParticleContainer.remap();

Sorting particles
-----------------

Sorting by cell
_______________

The particles can also be sorted by cell by calling ``void ParticleContainer<P>::sortParticleByCell()``, like so:

.. code-block:: c++

   ParticleContainer<P> myParticleContainer;

   myParticleContainer.sortParticlesByCell();

Internally in ``ParticleContainer<P>``, this will place the particles in another container which can be iterated over on a per-cell basis.
This is different from ``List<P>`` and ``ListBox<P>`` above, which contained particles stored on a per-patch basis with no internal ordering of the particles.

The per-cell particle container is a ``Vector<RefCountedPtr<LayoutData<BinFab<P> > > >`` type where again the ``Vector`` holds the particles on each AMR level and the ``LayoutData<BinFab>`` holds one ``BinFab`` on each grid patch.
The ``BinFab`` is also a template, and it holds a ``List<P>`` in each grid cell.
Thus, this data structure stores the particles per cell rather than per patch.
Due to the horrific template depth, this container is typedef'ed as ``AMRCellParticles<P>``.

To get cell-sorted particles one can call

.. code-block:: c++

   AMRCellParticles<P>& cellSortedParticles = myParticleContainer.getCellParticles();

Iteration over cell-sorted particles is mostly the same as for patch-sorted particles, except that we also need to explicitly iterate over the grid cells in each grid patch:

.. code-block:: c++

   const int comp = 0;

   // Iterate over all AMR levels
   for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){

      // Get the grids on this level
      const DisjointBoxLayout& dbl = m_amr->getGrids(myParticleContainer.getRealm())[lvl];

      // Iterate over grid patches on this level
      for (DataIterator dit(dbl); dit.ok(); ++dit){

         // Get the Cartesian box for the current grid aptch
         const Box cellBox = dbl[dit()];

	 // Get the particles in the current grid patch.
	 BinFab<P>& cellSortedBoxParticles = (*cellSortedParticles[lvl])[dit()];

	 // Iterate over all cells in the current box
	 for (BoxIterator bit(cellBox); bit.ok(); ++bit){
	    const IntVect iv = bit();

	    // Get the particles in the current grid cell.
	    List<P>& cellParticles = cellSortedBoxParticles(iv, comp);

	    // Do something with cellParticles
	    for (ListIterator<P> lit(cellParticles); lit.ok(); ++lit){
	       P& p = lit();
	    }
	 }
      }
   }

Sorting by patch
________________

If the particles need to return to patch-sorted particles:

.. code-block:: c++

   ParticleContainer<P> myParticleContainer;

   myParticleContainer.sortParticlesByPatch();

.. important::
   
   If particles are sorted by cell, calling ``ParticleContainer<P>`` member functions that fetch particles by patch will issue an error.
   This is done by design since the patch-sorted particles have been moved to a different container.
   Note that remapping particles also requires that the particles are patch-sorted.
   Calling ``remap()`` with cell-sorted particles will issue a run-time error. 

Allocating particles
--------------------

``AmrMesh`` has a very simple function for allocating a ``ParticleContainer<P>``:

.. code-block:: c++

  template <typename P>
  void allocate(ParticleContainer<P>& a_container, const int a_pvrBuffer, const std::string a_realm);		

which will allocate a ``ParticleContainer`` on realm ``a_realm`` with a buffer zone of ``a_pvrBuffer``. 
This buffer zone adjusts if particles on the fine side of a refinement boundary map to the coarse grid or the fine grid (see :ref:`Chap:ParticleMapping`). 

.. _Chap:ParticleMapping:   
   
Mapping and remapping
---------------------

Mapping particles with ParticleValidRegion
__________________________________________

The ``ParticleValidRegion`` (PVR) allows particles to be transferred to coarser grid levels if they are within a specified number of grid cells from the refinement boundary.
There are two reasons why such a functionality is useful:

#. Particles that live in the first strip of cells on the fine side of a refinement boundary have deposition clouds that hang over the boundary and into ghost cells.
   This mass must be added to the coarse level, which adds algorithmic complexity (``chombo-discharge`` can handle this complexity). 
   
#. Deposition and interpolation kernels can be entirely contained within a grid level.
   It might be useful to keep the kernel on a specific AMR level for a certain number of time step. 

.. figure:: /_static/figures/ParticleValidRegion.png
   :width: 50%
   :align: center

   The ``ParticleValidRegion`` allows particles whose position fall into a fine grid patch to be moved to a coarser level if they are within a specified distance from the refinement boundary.
   In this case, the green particles that overlap with the fine-level grid are remapped to the coarse level. 

The PVR is automatically allocated through the particle constructor by specifying the ``a_pvrBuffer`` flag.
If you do not want to use PVR functionality, simply set ``a_pvrBuffer = 0`` for your ``ParticleContainer<P>``.
In this case the particles will live on the grid patch that contains them. 


Remapping particles
___________________

Particles that move off their original grid patch must be remapped in order to ensure that they are assigned to the correct grid.
The remapping function for ``ParticleContainer<P>`` is ``void ParticleContainer<P>::remap()``, which is simply used as follows:
   
.. code-block::

   ParticleContainer<P> myParticles;

   myParticles.remap();

Note that if a PVR region is set, the particle container remapping will respect it. 

Regridding
----------

``ParticleContainer<P>`` is comparatively simple to regrid, and this is done in two steps:

1. Each MPI rank collects *all* particles on a single ``List<P>`` by calling

   .. code-block:: c++

      void ParticleContainer<P>::preRegrid(int a_base)
      
   This will pull the particles off their current grids and collect them in a single list (on a per-rank basis).
   
2. When ``ParticleContainer<P>`` regrids, each rank adds his ``List<P>`` back into the internal particle containers.

The use case typically looks like this:

.. code-block:: c++
   
   ParticleContainer<P> myParticleContainer;

   // Each rank caches his particles
   const int baseLevel = 0;
   myParticleContainer.preRegrid(0);

   // Driver does a regrid.
   .
   .
   .
   
   // After the regrid we fetch grids from AmrMesh:
   Vector<DisjointBoxLayout> grids;
   Vector<ProblemDomain> domains;
   Vector<Real> dx;
   Vector<int> refinement_ratios;
   int base;
   int newFinestLevel;
   
   myParticleContainer.regrid(grids, domains, dx, refinement_ratios, baseLevel, newFinestLevel);

Here, ``baseLevel`` is the finest level that didn't change and ``newFinestLevel`` is the finest AMR level after the regrid. 

.. _Chap:MaskedParticles:

Masked particles
----------------

``ParticleContainer<P>`` also supports the concept of *masked particles*, where one can fetch a subset of particles that live only in specified regions in space.
Typically, this "specified region" is the refinement boundary, but the functionality is generic and might prove useful also in other cases.

When *masked particles* are used, the user can provide a boolean mask over the AMR hierarchy and obtain the subset of particles that live in regions where the mask evaluates to true.
This functionality is for example used for some of the particle deposition methods in ``chombo-discharge`` where we deposit particles that live near the refinement boundary with special deposition functions.

To fill the masked particles, ``ParticleContainer<P>`` has members functions for copying the particles into internal data containers which the user can later fetch.
The function signatures for these are

.. code-block:: c++

   using AmrMask = Vector<RefCountedPtr<LevelData<BaseFab<bool> > > >;

   template <class P>
   void copyMaskParticles(const AmrMask& a_mask) const;

   template <class P>   
   void copyNonMaskParticles(const AmrMask& a_mask) const;

The argument ``a_mask`` holds a bool at each cell in the AMR hierarchy.
Particles that live in cells where ``a_mask`` is true will be copied to an internal data holder in ``ParticleContainer<P>`` which can be retried through a call

.. code-block:: c++

   AMRParticles<P>& maskParticles = myParticleContainer.getMaskParticles();

Note that ``copyNonMaskParticles`` is just like ``copyMaskParticles`` except that the bools in ``a_mask`` have been flipped.

Note that the mask particles are *copied*, and the original particles are left untouched.
After the user is done with the particles, they should be deleted through the functions ``void clearMaskParticles()`` and ``void clearNonMaskParticles``, like so:

.. code-block:: c++

   AmrMask myMask;
   ParticleContainer<P> myParticles;

   // Copy mask particles
   myParticles.copyMaskParticles(myMask);

   // Do something with the mask particles
   AMRParticles<P>& maskParticles = myParticleContainer.getMaskParticles();

   // Release the mask particles
   myParticles.clearMaskParticles();

Creating particle halo masks
____________________________

``AmrMesh`` can register a *halo* mask with a specified width:

.. code-block:: c++

   void registerMask(const std::string a_mask, const int a_buffer, const std::string a_realm);

where ``a_mask`` must be ``"s_particle_halo"``.
This will register a mask which is false everywhere except in coarse-grid cells that are within a distance a_buffer from the refinement boundary, see :numref:`Fig:HaloMask`.

.. _Fig:HaloMask:
.. figure:: /_static/figures/HaloMask.png
   :width: 40%
   :align: center

   Example of a particle halo mask (shaded green color) surrounding refined grid levels.
   
Embedded boundaries
-------------------

``ParticleContainer<P>`` is EB-agnostic and has no information about the embedded boundary.
This means that particles remap just as if the EB was not there.
Interaction with the EB is done via the implicit function or discrete information, as well as modifications in the interpolation and deposition steps. 

Signed distance function
________________________

When signed distance functions are used, one can always query how far a particle is from a boundary:

.. code-block:: c++

   List<P>& particles;
   BaseIF distanceFunction;

   for (ListIterator<P> lit(particles); lit.ok(); ++lit){
      const P& p          = lit();
      const RealVect& pos = p.position();

      const Real distanceToBoundary = distanceFunction.value(pos);
   }

If the particle is inside the EB then the signed distance function will be positive and the particle can be removed from the simulation.
The distance function can also be used to detect collisions between particles and the EB. 

Particle depositon
------------------

To deposit particles on the mesh, the user can call the templated function ``AmrMesh::depositParticles`` which has a signature

.. code-block:: c++
		
  template <class P, const Real&(P::*particleScalarField)() const>
  void depositParticles(EBAMRCellData&              a_meshData,
			const std::string&          a_realm,
			const phase::which_phase&   a_phase,	       
			const ParticleContainer<P>& a_particles,
			const DepositionType        a_depositionType,
			const CoarseFineDeposition  a_coarseFineDeposition,
			const bool                  a_forceIrregNGP);

  template <class P, const RealVect&(P::*particleVectorField)() const>
  void depositParticles(EBAMRCellData&              a_meshData,
			const std::string&          a_realm,
			const phase::which_phase&   a_phase,	       
			const ParticleContainer<P>& a_particles,
			const DepositionType        a_depositionType,
			const CoarseFineDeposition  a_coarseFineDeposition,
			const bool                  a_forceIrregNGP);			

Here, the template parameter ``P`` is the particle type and the template parameter ``particleScalarField`` is a C++ pointer-to-member-function.
This function must have the indicated signature ``const Real& P::particleScalarField() const`` *or* the signature ``Real P::particleScalarField() const``.
The pointer-to-member ``particleScalarField`` indicates the variable to be deposited on the mesh.
This function pointer does not need to return a member in the particle class.

Note that when depositing vector-quantities (such as electric currents), one must call the version which takes ``RealVect P::particleVectorField() const`` as a template parameter.
The supplied function must return a ``RealVect`` and ``a_meshData`` must have ``SpaceDim`` components. 

Next, the input arguments to ``depositParticles`` are the output mesh data holder (must have exactly one or ``SpaceDim`` components), the realm and phase where the particles live, and the particles themselves (``a_particles``).
The enum ``DepositionType`` and input argument ``a_depositionType`` indicates the deposition method.
Valid arguments are

* ``DepositionType::NGP`` (Nearest grid-point).
* ``DepositionType::CIC`` (Cloud-In-Cell).
* ``DepositionType::TSC`` (Triangle-Shaped Cloud).
* ``DepositionType::W4``  (Fourth order weighted).

The input argument ``a_coarseFineDeposition`` determines how coarse-fine deposition is handled.
Strictly speaking, this argument only affects how the particle mass is deposited from the coarse level to the fine level. 
Valid input arguments are

* ``CoarseFineDeposition::PVR`` This uses a standard PVR formulation.
  When the particles near the refinement boundary deposit on the mesh, some of the mass from the coarse-side particles will end up underneath the fine grid.
  This mass is interpolated to the fine grid using piecewise constant interpolation.
  If the fine-level particles also have particle clouds that hang over the refinement boundary, the hanging mass will be added to the coarse level.
* ``CoarseFineDeposition::Halo`` This uses a what we call *halo* particles. 
  Instead of interpolating the mass from the invalid coarse region onto the fine level, the particles near the refinement boundary (i.e., the *halo* particles) deposit directly into the fine level but with 2x or 4x the particle width.
  So, if a coarse-level particle lives right next to the fine grid and the refinement factor between the grids is :math:`r`, it will deposit both into the fine grid with :math:`r` times the particle width compared to the coarse grid.
  Again, if the fine-level particles also have particle clouds that hang over the refinement boundary, the hanging mass will be added to the coarse level.
* ``CoarseFineDeposition::HaloNGP`` This uses halo particles, but the particles along the refinement boundary are deposited with an NGP scheme. 

Finally, the flag ``a_forceIrregNGP`` permits the user to enforce nearest grid-point deposition in cut-cells.
This option is motivated by the fact that some applications might require hard mass conservation, and the user can ensure that mass is never deposited into covered grid cells. 

As an example, if the particle type ``P`` needs to deposit a computational mass on the mesh, the particle class will at least contain the following member functions:

.. code-block:: c++

   class P : public BinItem {
   public:

      const Real& mass() const {
         return m_mass;
      }

      Real mass2() const {
         return m_mass*m_mass.
      }

      RealVect momentum() const {
         return m_mass*m_velocity;
      }

   protected:

      Real m_mass;

      Real m_velocity;
   };

Here, we have included an extra member function ``mass()`` which returns the squared mass.
Note that the function does not return a member variable but an r-value.
When depositing the mass on the mesh the user will e.g. call

.. code-block:: c++

   RefCountedPtr<AmrMesh> amr;

   amr->depositParticles<P, &P::mass >(...);
   amr->depositParticles<P, &P::mass2>(...);

When depositing momentum, use

.. code-block:: c++
		
   amr->depositParticles<P,  &P::momentum>(...).      

Particle interpolation
----------------------

To interpolate a field onto a particle position, the user can call the ``AmrMesh`` member functions


.. code-block:: c++
		
  template <class P, Real&(P::*particleScalarField)()>
  void interpolateParticles(ParticleContainer<P>&      a_particles,
			    const std::string&         a_realm,
			    const phase::which_phase&  a_phase,	       			    
			    const EBAMRCellData&       a_meshScalarField,
			    const DepositionType       a_interpType,	
			    const bool                 a_forceIrregNGP) const;

  template <class P, RealVect&(P::*particleVectorField)()>
  void interpolateParticles(ParticleContainer<P>&      a_particles,
			    const std::string&         a_realm,
			    const phase::which_phase&  a_phase,	       			    			    
			    const EBAMRCellData&       a_meshVectorField,
			    const DepositionType       a_interpType,
			    const bool                 a_forceIrregNGP) const;

The function signature for particle interpolation is pretty much the same as for particle deposition, with the exception of the interpolated field.
The template parameter ``P`` still indicates the particle type, but the user can interpolate onto either a scalar particle variable or a vector variable.
For example, in order to interpolate the particle acceleration, the particle class (let's call it ``MyParticleClass``) will typically have a member function ``RealVect& acceleration()``, and in this case one can interpolate the acceleration by

.. code-block:: c++

   RefCountedPtr<AmrMesh> amr;

   amr->interpolateParticles<MyParticleClass, &MyParticleClass::acceleration>(...)

Note that if the user interpolates onto a scalar variable, the mesh variable must have exactly one component.
Likewise, if interpolating a vector variable, the mesh variable must have exact ``SpaceDim`` components.

.. _Chap:SuperParticles:

Superparticles
--------------

Custom approach
_______________

For a custom approach of managing superparticles, users can simply manipulate the particle lists in the grid patches or grid cells.
In each case one starts with a list ``List<P>`` that needs to be modified. 

kD-trees
________

Overview
^^^^^^^^

``chombo-discharge`` has functionality for spatially partitioning particles using kD-trees, which can be used as a basis for particle merging and splitting.
kD-trees operate by partitioning a set of input primitives into spatially coherent subsets.
At each level in the tree recursion one chooses an axis for partitioning one subset into two new subsets, and the recursion continues until the partitioning is complete.
:numref:`Fig:PartitionKD` shows an example where a set of initial particles are partitioned using such a tree. 

.. _Fig:PartitionKD:
.. figure:: /_static/figures/PartitionKD.png
   :width: 75%
   :align: center

   Example of a kD-tree partitioning of particles in a single cell.

.. note:: 

   The source code for the kD-tree functionality is given in :file:`$DISCHARGE_HOME/Source/Particle/CD_SuperParticles.H`.       

Particle partitioners
^^^^^^^^^^^^^^^^^^^^^

The kD-tree partitioner requires a user-supplied criterion for particle partitioning.
Only the partitioner ``PartitionEqualWeight`` is currently supported, and this partitioner will divide the original subset into two new subsets such that the particle weights in the two halves differs by at most one physical particle.
This partitioner is imlemented as

.. code-block:: c++

  template <class P, Real& (P::*weight)(), const RealVect& (P::*position)() const>
  typename KDNode<P>::Partitioner PartitionEqualWeight;

Here, ``P`` is the particle type, and this class *must* have function members ``Real& P::weight()`` and ``const RealVect& P::position()`` which return the particle weight and position.

.. warning::
   
   ``PartitionEqualWeight`` will usually split particles to ensure that the weight in the two subsets are the same (thus creating new particles). 
   In this case any other members in the particle type are copied over into the new particles.

The particles in each leaf of the kD-tree can then be merged into new particles.
Since the weight in the nodes of the tree differ by at most one, the resulting computational particles also have weights that differ by at most one.

.. _Fig:SuperKD:
.. figure:: /_static/figures/SuperKD.png
   :width: 50%
   :align: center

   kD-tree partitioning of particles into new particles whose weight differ by at most one.
   Left: Original particles with weights between 1 and 100.
   Right: Merged particles.
   

.. _Chap:ParticleOps:

ParticleOps
-----------

``ParticleOps`` is a static data class that provides methods for commonly used particle operations.
These include

* Intersection of particles with EBs.
* Intersection of particles with domain edges/faces.
* Drawing particles from a probability distribution.

The ``ParticleOps`` API is found at `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classParticleOps.html>`_.
