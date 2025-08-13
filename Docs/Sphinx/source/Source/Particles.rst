.. _Chap:Particles:

Particles
=========

``chombo-discharge`` supports computational particles using native ``Chombo`` particle data.
The source code for the particle functionality resides in :file:`$DISCHARGE_HOME/Source/Particle`.
Particle support contains the following basic features:

* Particle-mesh operations, i.e., deposition and interpolation of particle variables to/from the mesh.
* Particle distribution and remapping with MPI.
* Rudimentary particle output to H5Part files.

Particle support is generally speaking templated, so that users can define new particle types that contain a desired set of variables.
Typically, these will are derived from :ref:`Chap:GenericParticle`, which is discussed below.

.. _Chap:GenericParticle:

GenericParticle
---------------

``GenericParticle`` is a default particle usable by the ``Chombo`` particle library.
The particle type is essentially a template

.. code-block:: c++

   template <size_t M, size_t N>
   class GenericParticle
   {
   public
      RealVect&
      position();
   protected:
      RealVect m_position;
      std::array<Real, M> m_scalars;
      std::array<RealVcet, N> m_vectors;      
   };

where ``M`` and ``N`` are the number of ``Real`` and ``RealVect`` variables for the particle.
The ``GenericParticle`` always stores the position of the particle, which is available through ``GenericParticle<M,N>::position``.

To fetch the ``Real`` and ``RealVect`` variables, ``GenericParticle`` has member functions

.. code-block:: c++

   template <size_t K>
   inline Real&
   GenericParticle<M,N>::real();

   template <size_t K>
   inline RealVect&
   GenericParticle<M,N>::vect();

If using ``GenericParticle`` directly, the correct C++ way of fetching one of these variables is

.. code-block:: c++
		
   GenericParticle<2,2> p;

   Real& s = p.template real<0>();

Note that one must include the template keyword.    

.. tip::

   The ``GenericParticle`` C++ API is found at `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classGenericParticle.html>`_.

Custom particles
----------------

To create a simple custom particle class with more sane signatures, one can inherit from ``GenericParticle`` and specify new function signatures that return the appropriate fields.
An example of this is given in the code-block below, where we define ``KineticParticle`` to be a particle that contains the th5ree additional fields on top of ``GenericParticle`` (weight, velocity, and acceleration).

.. code-block:: c++
	   
   class KineticParticle : public GenericParticle<1,2>
   {
   public:
      inline
      Real& weight() {
         return this->real<0>();
      }
      
      inline
      RealVect& velocity() {
         return this->vect<0>();
      }

      inline
      RealVect& acceleration() {
         return this->vect<1>();
      }            
   };

There are many particles in ``chombo-discharge``, see the ``GenericParticle`` C++ API for more information. 

.. _Chap:ParticleContainer:

ParticleContainer
------------------

The ``ParticleContainer<P>`` is a template class that

#. Stores computational particles of type ``P`` over an AMR hierchy.
#. Provides infrastructure for remapping particles.
#. Provides functionality for getting a list of particles within a specified grid patch.
#. Provides functionality that is required during regrids.
#. Other types of functionality, like grouping particles into grid cells, set and get functions for assigning particle variables, etcl.

``ParticleContainer<P>`` uses the ``Chombo`` structure ``ParticleData<P>`` under the hood, and therefore has template constraints on ``P``.
The simplest way to use ``ParticleContainer`` for a new type of particle is to let ``P`` inherit from :ref:`Chap:GenericParticle`, which will fulfill all template constraints.

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
``ParticleData<P>`` can be thought of as a ``LevelData<ListBox<P>>``, although it actually inherits from ``LayoutData<ListBox<P>>``.
Each grid patch contains a ``ListBox<P>`` of particles.


AMRParticles<P>
_______________

``AMRParticles<P>`` is our AMR version of ``ParticleData<P>``.
It is a simply a typedef of a vector of pointers to ``ParticleData<P>`` on each level:

.. code-block:: c++

   template <class P>
   using AMRParticles = Vector<RefCountedPtr<ParticleData<P>>>;

Again, the ``Vector`` indicates the AMR level and the ``ParticleData<P>`` is a distributed data holder that holds the particles on each AMR level.

``AMRParticles<P>`` always lives within ``ParticleContainer<P>``, and is the class member of ``ParticleContainer<P>`` that actually holds the particles.

Basic usage
-----------

Here, we give some examples of basic usage of ``ParticleContainer``.
For the full API, see the ``ParticleContainer`` C++ API `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classParticleContainer.html>`_.

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

The per-cell particle container is a ``Vector<RefCountedPtr<LayoutData<BinFab<P>>>>`` type where again the ``Vector`` holds the particles on each AMR level and the ``LayoutData<BinFab>`` holds one ``BinFab`` on each grid patch.
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
  void allocate(ParticleContainer<P>& a_container, const std::string a_realm);		

which will allocate a ``ParticleContainer`` on realm ``a_realm``.
See :ref:`Chap:AmrMesh` for further details. 

.. _Chap:ParticleMapping:   
   
Particle mapping
----------------

Particles that move off their original grid patch must be remapped in order to ensure that they are assigned to the correct grid.
The remapping function for ``ParticleContainer<P>`` is

.. code-block:: c++

   void
   ParticleContainer<P>::remap();

This is simply used as follows:
   
.. code-block::

   ParticleContainer<P> myParticles;

   myParticles.remap();

Regridding
----------

``ParticleContainer<P>`` is comparatively simple to regrid, and this is done in two steps:

1. *Before* creating the new grids, each MPI rank collects *all* particles on a single ``List<P>`` by calling

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

   using AmrMask = Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>;

   template <class P>
   void copyMaskParticles(const AmrMask& a_mask) const;

   template <class P>   
   void copyNonMaskParticles(const AmrMask& a_mask) const;

The argument ``a_mask`` holds a bool at each cell in the AMR hierarchy.
Particles that live in cells where ``a_mask`` is true will be copied to an internal data holder in ``ParticleContainer<P>`` which can be retrieved through a call

.. code-block:: c++

   AMRParticles<P>& maskParticles = myParticleContainer.getMaskParticles();

Note that ``copyNonMaskParticles`` is just like ``copyMaskParticles`` except that the bools in ``a_mask`` have been flipped.

In the above functions the mask particles are *copied*, and the original particles are left untouched.
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
This will register a mask which is false everywhere except in coarse-grid cells that are within a distance ``a_buffer`` from the refinement boundary, see :numref:`Fig:HaloMask`.
This functionality is useful when processing particles on the refinement boundary using special deposition functions since the halo mask allows us to straightforwardly extract those particles.

.. _Fig:HaloMask:
.. figure:: /_static/figures/HaloMask.png
   :width: 40%
   :align: center

   Example of a particle halo mask (shaded green color) surrounding refined grid levels.

.. _Chap:ParticleEB:
   
Wall interaction
----------------

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
See :ref:`Chap:AmrMesh` for details on how to obtain the distance function.

Domain edges
____________

By default, the ``ParticleContainer`` remapping function will discard particles that fall outside of the domain.
The user can also check if this happen by checking if the particle position is outside the computational domain:

.. code-block:: c++

   GenericParticle<0,0> p;

   const RealVect pos    = p.position();
   const RealVect probLo = m_amr->getProbLo();
   const RealVect probHi = m_amr->getProbHi();

   bool outside = false;
   for (int dir = 0; dir < SpaceDim; dir++) {
      if(pos[dir] < probLo[dir] || pos[dir] > probHi[dir]) {
         outside = true;
      }
   }

Particle intersection
_____________________

It is occasionally useful to catch particles that hit an EB or crossed a domain side.
Assuming that the particle type ``P`` also has a member function that stores the starting position of the particle, one can compute the intersection point between the particle trajectory and the EB and domain edges/faces.
Currently, :ref:`Chap:AmrMesh` supports two methods for computing this

* Using a bisection algorithm with a user-specified step.
* Using a ray-casting algorithm.

These algorithms differ in the sense that the bisection approach will check for a particle crossing between two positions :math:`\mathbf{x}_0` and :math:`\mathbf{x}_1` using a pre-defined tolerance.
The ray-casting algorithm will check if the particle can move from :math:`\mathbf{x}_0` towards :math:`\mathbf{x}_1` by using a variable step along the particle trajectory.
This step is selected from the signed distance from the particle position to the EB such that it uses a large step if the particle is far away from the EB.
Conversely, if the particle is close to the EB a small step will be used.
For the function signatures, see :ref:`Chap:AmrParticleIntersection`.
The algorithms that operate under the hood of these routines are given in ``ParticleOps``, see :ref:`Chap:ParticleOps`. 

Both the bisection and ray-casting algorithm have weaknesses.
The bisection algorithm algorithm requires a user-supplied step in order to operate efficiently, while the ray-casting algorithm is very slow when the particle is close to the EB and moves tangentially along it.
Future versions of ``chombo-discharge`` will likely include more sophisticated algorithms. 

.. tip::

   ``AmrMesh`` also stores the implicit function on the mesh, which could also be used to resolved particle collisions with the EB/domain.

.. _Chap:ParticleMesh:

Particle-mesh
-------------

Particle deposition
___________________

To deposit particles on the mesh, the user can call the templated function ``AmrMesh::depositParticles`` which have a signatures

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

When depositing vector-quantities (such as electric currents), one must call the version which takes ``RealVect P::particleVectorField() const`` as a template parameter.
The supplied function must return a ``RealVect`` and ``a_meshData`` must then have ``SpaceDim`` components. 

Next, the input arguments to ``depositParticles`` are the output mesh data holder (must have exactly one or ``SpaceDim`` components), the realm and phase where the particles live, and the particles themselves (``a_particles``).
Finally, the flag ``a_forceIrregNGP`` permits the user to enforce nearest grid-point deposition in cut-cells.
This option is motivated by the fact that some applications might require hard mass conservation, and the user can ensure that mass is never deposited into covered grid cells.

The input argument ``a_depositionType`` indicates the deposition method, while ``a_coarseFineDeposition`` deposition modifications near refinement boundaries.
These are discussed below.

Base deposition
^^^^^^^^^^^^^^^

The base deposition scheme is specified by an enum ``DepositionType`` with valid values:

* ``DepositionType::NGP`` (Nearest grid-point).
* ``DepositionType::CIC`` (Cloud-In-Cell).
* ``DepositionType::TSC`` (Triangle-Shaped Cloud).
* ``DepositionType::W4``  (Fourth order weighted).


Coarse-fine deposition
^^^^^^^^^^^^^^^^^^^^^^

The input argument ``a_coarseFineDeposition`` determines how coarse-fine deposition is handled.
Refinement boundaries introduce additional complications in the deposition scheme due to

#. Fine-grid particles whose deposition clouds hang over the refinement boundary and onto the coarse level.
#. Coarse-grid particles whose deposition clouds stick underneath the fine-level.

.. _Fig:ParticleDeposition:
.. figure:: /_static/figures/ParticleDeposition.png
   :width: 75%
   :align: center

   Sketch of deposition schemes near refinement boundaries and cut-cells. 

``chombo-discharge`` currently supports three methods for handling coarse-fine deposition.
In all of these methods the mass on the fine grid particles whose deposition clouds hang over the refinement boundaries is simply added to the coarse grid.
For the coarse-grid particles the following processes then occur:

* ``CoarseFineDeposition::Interp`` 
  This method permits the coarse-grid particles to deposit into the region underneath the fine grid.
  The deposited quantity is then piecewise interpolated onto the fine grid.
  The indicated coarse-grid deposition cloud in :numref:`Fig:ParticleDeposition` will then add its mass into two layers of fine-grid cells. 
* ``CoarseFineDeposition::Halo``
  This method extracts the coarse-grid particles that live on the refinement boundary.
  These particles are then transferred to the fine level and they are then deposit on the fine grid using the original particle width.
  For example, if using a CIC scheme and having a refinement factor of 2 between the coarse grid and fine grid, the particle width on the fine grid will be :math:`2\Delta x_{\textrm{fine}}` rather then :math:`\Delta x_{\textrm{fine}}`.

  .. warning::

     This functionality is currently limited to NGP and CIC schemes.
     
* ``CoarseFineDeposition::HaloNGP``
  Like ``CoarseFineDeposition::Halo``, this method also extracts the coarse-grid particles on the coarse side of the refinement boundary, but rather than using the original scheme these particles are deposited with an NGP scheme. 



Particle interpolation
______________________

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

Example
_______

Assume that we have some particle class ``KineticParticle`` defined as

.. code-block:: c++
	   
   class KineticParticle : public GenericParticle<1,2>
   {
   public:
      inline
      Real& weight() {
         return this->real<0>();
      }
      
      inline
      RealVect& velocity() {
         return this->vect<0>();
      }

      inline
      RealVect& acceleration() {
         return this->vect<1>();
      }

      inline
      RealVect momentum() const {
         return this->weight() * this->velocity();
      }
   };

To deposit the weight, velocity, and momentum on the grid we would call

.. code-block:: c++

   RefCountedPtr<AmrMesh> amr;

   amr->depositParticles<KineticParticle, &KineticParticle::mass >(...);
   amr->depositParticles<KineticParticle, &KineticParticle::velocity>(...);
   amr->depositParticles<KineticParticle, &KineticParticle::momentum>(...);

Likewise, to interpolate onto these fields we can call

.. code-block:: c++

   RefCountedPtr<AmrMesh> amr;

   amr->interpolateParticles<KineticParticle, &KineticParticle::mass >(...);
   amr->interpolateParticles<KineticParticle, &KineticParticle::velocity>(...);

.. _Chap:ParticleVisualization:   

Particle visualization
----------------------

.. note::

   Particle visualization is currently a work in progress.

Simple particle visualization can be performed by writing ``H5Part`` compatible files which can be read by VisIt.
This is done through the function ``writeH5Part`` in the ``DischargeIO`` namespace, with the following signature:

.. code-block:: c++

  template <size_t M, size_t N>
  void
  writeH5Part(const std::string                               a_filename,
              const ParticleContainer<GenericParticle<M, N>>& a_particles,
              const std::vector<std::string>                  a_realVars = std::vector<std::string>(),
              const std::vector<std::string>                  a_vectVars = std::vector<std::string>(),
              const RealVect                                  a_shift    = RealVect::Zero) noexcept;

This routine permits particles to be written (in parallel, when using MPI) into a file readable by VisIt.
While users will typically not work directly with ``GenericParticle``, casting to a proper format is quite simple, e.g.

.. code-block:: c++

   // ItoParticle inherits GenericParticle<5,3>		
   ParticleContainer<ItoParticle> myParticles;

   DischargeIO::writeH5Part("my_particles.h5part", (const ParticleContainer<GenericParticle<5,3>>&) myParticles);

The optional arguments ``a_realVars`` and ``a_vectVars`` permit the user to set the output variable names for the ``M`` scalar variables and the ``N`` vector variables.
The argument ``a_shift`` will simply shift the particle positions in the output HDF5 file. 

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
   :width: 90%
   :align: center

   Example of a kD-tree partitioning of particles in a single cell.

.. tip:: 

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
   :width: 75%
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
