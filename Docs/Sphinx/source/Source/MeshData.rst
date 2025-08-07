.. _Chap:MeshData:

Mesh data
=========

Mesh data structures of the type discussed in :ref:`Chap:SpatialDiscretization` are derived from a class ``EBAMRData<T>`` which holds a ``T`` in every grid patch across the AMR hiearchy.
A requirement on the datatype ``T`` is that it must be linearizable so that it can be communicated across MPI ranks. 
Internally, the data is stored as a ``Vector<RefCountedPtr<LevelData<T>>>``.
Here, the ``Vector`` holds data on each AMR level; the data is allocated with a smart pointer called ``RefCountedPtr`` which points to a ``LevelData`` template structure, see :ref:`Chap:Basics`.
The first entry in the ``Vector`` is base AMR level and finer levels follow later in the ``Vector``.
:numref:`Fig:EBAMRData` shows a sketch of the data layout in a two-level AMR hierarchy. 

.. _Fig:EBAMRData:
.. figure:: /_static/figures/PatchBasedAMR.png
   :width: 60%
   :align: center

   Cartesian patch-based refinement showing two grid levels.
   The finer level consists of two patches (red and blue zones), and these zones each have a ghost cell layer 2 cells wide.
   The data lies on top of a coarse-grid data, i.e., data simultaneously exists on both the fine and the coarse levels. 
   This data type is encapsulated by ``EBAMRData<T>``.

The reason for having class encapsulation of mesh data is due to :ref:`Chap:Realm`, so that we can only keep track on which ``Realm`` the mesh data is defined.
Users will interact with ``EBAMRData<T>`` through application code, or interacting with the core AMR functionality in :ref:`Chap:AmrMesh` (such as computing gradients, interpolating ghost cells etc.).
:ref:`Chap:AmrMesh` has functionality for constructing most ``EBAMRData<T>`` types on a ``Realm``, and ``EBAMRData<T>`` itself it typically not used anywhere elsewhere within ``chombo-discharge``.

A number of explicit template specifications exist and are frequently used.
These are outlined below: 

.. code-block:: c++

   typedef EBAMRData<EBCellFAB>        EBAMRCellData;  // Cell-centered single-phase data
   typedef EBAMRData<EBFluxFAB>        EBAMRFluxData;  // Face-centered data in all coordinate direction
   typedef EBAMRData<EBFaceFAB>        EBAMRFaceData;  // Face-centered in a single coordinate direction
   typedef EBAMRData<BaseIVFAB<Real> > EBAMRIVData;    // Data on irregular data centroids
   typedef EBAMRData<DomainFluxIFFAB>  EBAMRIFData;    // Data on domain phases
   typedef EBAMRData<BaseFab<bool> >   EBAMRBool;      // For holding bool at every cell

   typedef EBAMRData<MFCellFAB>        MFAMRCellData;  // Cell-centered multifluid data
   typedef EBAMRData<MFFluxFAB>        MFAMRFluxData;  // Face-centered multifluid data
   typedef EBAMRData<MFBaseIVFAB>      MFAMRIVData;    // Irregular face multifluid data   


For example, ``EBAMRCellData`` is a ``Vector<RefCountedPtr<LevelData<EBCellFAB> > >``, describing cell-centered data across the entire AMR hierarchy.
There are many more data structures in place, but the above data structures are the most commonly used ones.
Here, ``EBAMRFluxData`` is precisely like ``EBAMRCellData``, except that the data is stored on *cell faces* rather than cell centers.
Likewise, ``EBAMRIVData`` is a data holder that holds data on each EB centroid (or boundary centroid) across the entire AMR hierachy.
In the same way, ``EBAMRIFData`` holds data on each face of all cut-cells in the hierarchy. 

Allocating mesh data
--------------------

To allocate data over a particular ``Realm``, the user will interact with :ref:`Chap:AmrMesh`:

.. code-block:: c++

   const int numComps = 1;
   EBAMRCellData myData;
   m_amr->allocate(myData, "myRealm", phase::gas, numComps);

Here, ``numCOmps`` determine the number of data components.
Note that it *does* matter on which ``Realm`` and on which ``phase`` the data is defined.
See :ref:`Chap:Realm` for details.

The user *can* specify a number of ghost cells for his/hers application code directly in the ``AmrMesh::allocate`` routine, like so:

.. code-block:: c++

   int numComps = 1;
   EBAMRCellData myData;
   m_amr->allocate(myData, "myRealm", phase::gas, numComps, 5*IntVect::Unit);

If the user does not specify the number of ghost cells when calling ``AmrMesh::allocate``, :ref:`Chap:AmrMesh` will use the default number of ghost cells specified in the input file.

.. _Chap:MeshIteration:

Iterating over the AMR hierarchy
--------------------------------

To iterate over data in an AMR hierarchy, you will first iterate over levels and then the patches on each level:

.. code-block:: c++

   for (int lvl = 0; lvl < myData.size(); lvl++){
      LevelData<EBCellFAB>& levelData = *myData[lvl];

      const DisjointBoxLayout& levelGrids = levelData.disjointBoxLayout();
      
      for (DataIterator dit = levelGrids.dataIterator(); dit.ok(); ++dit){
         EBCellFAB& patchData = levelData[dit()];
      }
   }

Throughout ``chombo-discharge`` it will be common to see the above implemented explicitly as a loop that supports OpenMP:

.. code-block::

   for (int lvl = 0; lvl < myData.size(); lvl++){
      LevelData<EBCellFAB>& levelData = *myData[lvl];

      const DisjointBoxLayout& levelGrids = levelData.disjointBoxLayout();
      const DataIterator& dataIterator    = levelGrids.dataIterator();

      const int numBoxes = dataIterator.size();

   #pragma omp parallel for schedule(runtime)
      for (int currentBox = 0; currentBox < numBoxes; currentBox++) {
         const DataIndex& dataIndex = dataIterator[currentBox];
	 
         EBCellFAB& patchData = levelData[dataIndex];
      }
   }

Iterating over cells
--------------------

For single-valued data, ``chombo-discharge`` uses standard loops (in column-major order) for iterating over data.
For example, the standard loops for iterating over cell-centered data are

.. code-block:: c++

   namespace BoxLoops {
   
      template <typename Functor>
      ALWAYS_INLINE void
      loop(const Box& a_computeBox, Functor&& kernel, const IntVect& a_stride = IntVect::Unit);

      template <typename Functor>
      ALWAYS_INLINE void
      loop(VoFIterator& a_iter, Functor&& a_kernel);
   }

Here, the ``Functor`` argument is a C++ lambda or ``std::function`` which takes a grid cell as a single argument.
For the first loop, we iterate over all grid cells in ``a_computeBox``.
Iterating over the cut-cells in a patch data holder (like the ``EBCellFAB``) can be done with a ``VoFIterator``, which can iterate through cells on an ``EBCellFAB`` that are not covered by the geometry.
For example:

.. code-block:: c++

   const int component = 0;

   for (int lvl = 0; lvl < myData.size(); lvl++){
      LevelData<EBCellFAB>& levelData = *myData[lvl];

      const DisjointBoxLayout& levelGrids = levelData.disjointBoxLayout();
      
      for (DataIterator dit = levelGrids.dataIterator(); dit.ok(); ++dit){
         EBCellFAB& patchData       = levelData[dit()];
	 BaseFab<Real>& regularData = patchData.getSingleValuedFab();

	 auto regularKernel = [&](const IntVect& iv) -> void {
	    regularData(iv, component) = 1.0;
	 };

	 auto irregularKernel = [&](const VolIndex& vof) -> void {
	    patchData(vof, component = 1.0;
	 };

	 // Kernel regions (defined by user)
	 Box computeBox = ...
	 VoFIterator vofit = ...

	 BoxLoops::loop(computeBox, regularKernel);
	 BoxLoops::loop(vofit, irregularKernel);	 
      }
   }

There are loops available for other types of data (e.g., face-centered data), see the `BoxLoop documentation <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/CD__BoxLoops_8H.html>`_.



.. _Chap:Coarsening:

Coarsening data
---------------

Conservative coarsening of data is done using the ``averageDown(...)`` functions in :ref:`Chap:AmrMesh`.
When using these functions, coarse-grid data is replaced by a conservative average of fine grid data throughout the entire AMR hierarchy.
The signatures for various types of data are as follows:

.. code-block:: c++

   // Conservatively coarsen multifluid cell-centered data
   void averageDown(MFAMRCellData& a_data, const std::string a_realm) const;

   // Conservatively coarsen multifluid face-centered data
   void averageDown(MFAMRFluxData& a_data, const std::string a_realm) const;

   // Conservatively coarsen cell-centered data
   void averageDown(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase) const;

   // Conservatively coarsen face-centered data   
   void averageDown(EBAMRFluxData& a_data, const std::string a_realm, const phase::which_phase a_phase) const;

   // Conservatively coarsen EB-centered data      
   void averageDown(EBAMRIVData& a_data, const std::string a_realm, const phase::which_phase a_phase) const;  

There are other types of coarsening available also.
For example, the ``averageFaces(...)`` will use unweighted averaging, see the `AmrMesh API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classAmrMesh.html>`_ for further details. 

.. _Chap:GhostCells:

Filling ghost cells
-------------------

Filling ghost cells is done using the ``interpGhost(...)`` functions in :ref:`Chap:AmrMesh`.

.. code-block:: c++

   void interpGhost(MFAMRCellData& a_data, const std::string a_realm) const;

   void interpGhost(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase) const;

This will fill the specified number of ghost cells using data from the coarse level only, using piecewise linear interpolation. 

As an alternative, one *can* interpolate a single layer of ghost cells using the multigrid interpolator (see :ref:`Chap:MultigridInterpolation`).
In this case only a single layer of ghost cells are filled in regular regions, but additional ghost cells (up to some specified range) are filled near the EB.
This is often required when computing gradients (to avoid reaching into invalid cut-cells), see :ref:`Chap:Gradients` for details.
The functions for filling ghost cells in this way are

.. code-block:: c++

   void interpGhostMG(MFAMRCellData& a_data, const std::string a_realm) const;

   void interpGhostMG(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase) const;

See the `AmrMesh API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classAmrMesh.html>`_ for further details. 

.. _Chap:Gradients:

Computing gradients
-------------------

In ``chombo-discharge`` gradients are computed using a standard second-order stencil based on finite differences.
This is true everywhere except near the refinement boundary and EB where the coarse-side stencil will avoid using the coarsened data beneath the fine level.
This is shown in :numref:`Fig:EBGradient` which shows the typical 5-point stencil in regular grid regions, and also a much larger and more complex stencil.

In :numref:`Fig:EBGradient` we have shown two regular 5-point stencils (red and green).
The coarse stencil (red) reaches underneath the fine level and uses the data defined by coarsening of the fine-level data.
The coarsened data in this case is just an average of the fine-level data.
Likewise, the green stencil reaches over the refinement boundary and into one of the ghost cells on the coarse level.

:numref:`Fig:EBGradient` also shows a much larger stencil (blue stencil).
The larger stencil is necessary because computing the :math:`y` component of the gradient using a regular 5-point stencil would have the stencil reach underneath the fine level and into coarse data that is also irregular data.
Since there is no unique way (that we know of) for coarsening the cut-cell fine-level data onto the coarse cut-cell without introducing spurious artifacts into the gradient, we reconstruct the gradient using a least squares procedure.
In this case we fetch a sufficiently large neighborhood of cells for computing a least squares minimization of a local solution reconstruction in the neighborhood of the coarse cell.
In order to avoid fetching potentially badly coarsened data, this neighborhood of cells only uses *valid* grid cells, i.e. the stencil does not reach underneath the fine level at all.
Once this neighborhood of cells is obtained, we compute the gradient using the procedure in :ref:`Chap:LeastSquares`. 

.. _Fig:EBGradient:
.. figure:: /_static/figures/EBGradient.png
   :width: 50%
   :align: center

   Example of stencils for computing gradients near embedded boundaries.
   The red stencil shows a regular 5-point stencil for computing the gradient on the coarse side of the refinement boundary; it reaches into the coarsened data beneath the fine level.
   The green stencil shows a similar 5-point stencil on the fine side of the refinement boundary; the stencil reaches over the refinement boundary and into one ghost cell.
   The blue stencils shows a much more complex stencil which is computed using a least squares reconstruction procedure. 

To compute gradients of a scalar, one can simply call ``AmrMesh::computeGradient(...)`` functions:

.. code-block:: c++

  void computeGradient(EBAMRCellData&           a_gradient,
		       const EBAMRCellData&     a_phi,
                       const std::string        a_realm,
                       const phase::which_phase a_phase) const;

  void computeGradient(MFAMRCellData& a_gradient, const MFAMRCellData& a_phi, const std::string a_realm) const;		

See :ref:`Chap:AmrMesh` or refer to the `AmrMesh API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classAmrMesh.html>`_ for further details.

.. _Chap:CopyingData:

Copying data
------------

To copy data, one may use the ``EBAMRData<T>::copy(...)`` function *or* ``DataOps::copy`` (see :ref:`Chap:DataOps`).
These differ in the following way:

* ``EBAMRData<T>::copy`` works across realms, but will not copy ghost cells. 
* ``DataOps::copy`` will always do a local copy, and thus the data that is copied *must* be defined on the same realm.
  
If you call ``EBAMRData<T>::copy(...)``, the data holders will first check if they are both defined on the same realm.
If they are, a purely local copy is perform, which will include ghost cells. 
Communication copies involving MPI are performed otherwise, in which case ghost cells are *not* copied into the new data holder. 

.. _Chap:DataOps:

DataOps
-------

We have prototyped functions for many common data operations in a static class ``DataOps``.
For example, setting the value of various data holders can be done with

.. code-block:: c++

   EBAMRFluxData cellData;
   EBAMRFluxData fluxData;
   EBAMRIVData   irreData;
   
   DataOps::setValue(cellData, 0.0);
   DataOps::setValue(fluxData, 1.0);
   DataOps::setValue(irreData, 2.0);

For the full API, see the `DataOps documentation <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classDataOps.html>`_.   
