.. _Chap:ParticleData:

Understanding particle data
===========================

ParticleData
------------

On each grid level, particles live in templated data holders

.. code-block:: c++

   template <class T>
   ParticleData<T>

where ``T`` is the particle type.
Although ``ParticleData<T>`` can be thought of as a ``LevelData<T>``, it actually inherits ``LayoutData<ListBox<T> >``.
However, ``ParticleData<T>`` has routines for communicating particles with MPI.

On each patch, the particles are stored in a ``ListBox<T>``.
This class essentially consists of the domain box for the patch which holds the particles, and the particles themselves which are stored in a list ``List<T>``.
The particles are retrieved from a ``ListBox<T>`` object with a ``listItems()`` member function.
Here is an example:

.. code-block:: c++

   // Assume PD is our particle data and DBL is our grids
   for (DataIterator dit = DBL.dataIterator(); dit.ok(); ++dit){
      ListBox<T>& lb = PD[dit()];
      List<T>& particles = lb.listItems();

      // ... do something with the particles
   }

There are routines for allocating a particle data holder in ``amr_mesh``, which simply looks like this:

.. code-block:: c++
		
   template <typename T >
   void allocate(Vector<RefCountedPtr<ParticleData<T> > >& a_particles);

This will simply create the ``ParticleData<T>`` object on all the AMR levels (without any particles). 

Iterating over particles
________________________

In order to iterate over particles, you will use an iterator without random access:

.. code-block:: c++

   // Assume PD is our particle data and DBL is our grids
   for (DataIterator dit = DBL.dataIterator(); dit.ok(); ++dit){
      ListBox<T>& lb = PD[dit()];
      List<T>& particles = lb.listItems();

      ListIterator<T> lit(particles);
      for (lit.rewind(); lit; ++lit){
         T& p = particles[lit];

	 // ... do something with this particle
      }
   }


ParticleValidRegion
-------------------

The ``ParticleValidRegion`` (PVR) allows particles to be transferred to coarser grid levels if they are within a specified number of grid cells from the refinement boundary.
A compelling reason to do this is that if the particle lives on the refinement boundary, its deposition cloud will hang over the refinement boundary and into the ghost cells.
So, it is useful to keep the particles on the grid in such a way that the deposition and interpolation kernels are entirel contained within the grid.

.. figure:: figures/pvr.png
   :width: 480px
   :align: center

   The ParticleValidRegion allows particles whose position fall into a fine grid patch to be moved to a coarser level if they are within a specified distance from the refinement boundary. In this case, the green particles that overlap with the fine-level grid are placed in a ``ParticleData<T>`` holder on the coarse grid level.

To allocate a PVR,
Allocation of a PVR is done from ``amr_mesh`` (alternatively, call the constructor yourself) as follows:

.. code-block:: c++

   void allocate(EBAMRPVR& a_pvr, const int a_buffer); // buffer is the number of cells from the refinement boundary

Here, ``EBAMRPVR`` is simply a typedef'ed ``Vector<RefCountedPtr<ParticleValidRegion> >``. 

Remapping particles
-------------------

Particle remapping to the correct MPI ranks must be done if a particle leaves a grid patch and enters a different one, or leaves over a refinement boundary.
The figure below shows some typical cases.


.. figure:: figures/outcast.png
   :width: 480px
   :align: center

   Three cases of particle remapping. Here, the green particle stays on the fine level but needs to change MPI ownership.
   The red particle moves from the coarse level and to the fine level and needs to change both ownership and will also be deposited on the fine AMR level.
   The blue particle has moved to a different patch on the fine level but falls outside the fine level PVR must therefore be transferred to the coarse level.

Generally, remapping is done by first remapping particles on the level, and then transferring particles to coarser or finer grid level depending on whether or not they fall inside or outside the PVR.    

Outcasts
________

After particles have moved, ``ParticleData<T>`` has a method for locally gathering particles that are no longer in the correct grid patch, and another method for distributed the particles to the correct grid patch.
This is done as follows

.. code-block:: c++

   // Assume PD is a ParticleData<T> object

   PD.gatherOutcast();
   PD.remapOutcast();

We remark that

1. Some particles may remain in the outcast list
2. The remapping does not respect the PVR on each level.

Two-level transfer
__________________



Limitations
-----------
