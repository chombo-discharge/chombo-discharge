.. _Chap:MeshData:

Mesh data
=========

Mesh data structures in ``chombo-discharge`` are derived from a class ``EBAMRData<T>`` which holds a ``T`` in every grid patch across the AMR hiearchy.
Internally, the data is stored as a ``Vector<RefCountedPtr<LevelData<T> > >``.
Here, the ``Vector`` holds data on each AMR level; the data is allocated with a smart pointer called ``RefCountedPtr`` which points to a ``LevelData`` template structure, see :ref:`Chap:Basics`.
The first entry in the Vector is base AMR level and finer levels follow later in the Vector.

The reason for having class encapsulation of mesh data is due to :ref:`Chap:Realm`, so that we can only keep track on which ``Realm`` the mesh data is defined.
Users will not have to interact with ``EBAMRData<T>`` directly, but will do so primarily through application code, or interacting with the core AMR functionality in :ref:`Chap:AmrMesh` (such as computing gradients, interpolating ghost cells etc.).
``AmrMesh`` (see :ref:`Chap:AmrMesh`) has functionality for defining most ``EBAMRData<T>`` types one a ``Realm``, and ``EBAMRData<T>`` itself it typically not used anywhere elsewhere within ``chombo-discharge``.

A number of explicit template specifications also exist, and these are outlined below: 

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
Likewise, ``EBAMRIVData`` is a typedef'ed data holder that holds data on each cut-cell center across the entire AMR hierachy.
In the same way, ``EBAMRIFData`` holds data on each face of all cut cells. 

Allocating mesh data
--------------------

To allocate data over a particular ``Realm``, the user will interact with ``AmrMesh`` (see :ref:`Chap:AmrMesh`):

.. code-block:: c++

   int nComps = 1;
   EBAMRCellData myData;
   m_amr->allocate(myData, "myRealm", phase::gas, nComps);

Note that it *does* matter on which ``Realm`` and on which ``phase`` the data is defined.
See :ref:`Chap:Realm` for details.

The user *can* specify a number of ghost cells for his/hers application code directly in the ``AmrMesh::allocate`` routine, like so:

.. code-block:: c++

   int nComps = 1;
   EBAMRCellData myData;
   m_amr->allocate(myData, "myRealm", phase::gas, nComps, 5*IntVect::Unit);

If the user does not specify the number of ghost cells when calling ``AmrMesh::allocate``, ``AmrMesh`` will use the default number of ghost cells specified in the input file.   

Iterating over data
-------------------

To iterate over data in an AMR hierarchy, you will first iterate over levels and the patches in levels:

.. code-block:: c++

   for (int lvl = 0; lvl < myData.size(); lvl++){
      LevelData<EBCellFAB>& levelData = *myData[lvl];

      const DisjointBoxLayout& levelGrids = levelData.disjointBoxLayout();
      
      for (DataIterator dit = levelGrids.dataIterator(); dit.ok(); ++dit){
         EBCellFAB& patchData = levelData[dit()];
      }
   }



Iterating over the cells in a patch data holder (like the ``EBCellFAB``) can be done with a ``VoFIterator``.
The ``VoFIterator`` can iterate through cells on an ``EBCellFAB`` that are not covered by the geometry, but it can be slow to define so the typical iteration structure consists of a loop that iterates through all cells first, and then the irregular cells later.
For example:

.. code-block:: c++

   const int component = 0;

   for (int lvl = 0; lvl < myData.size(); lvl++){
      LevelData<EBCellFAB>& levelData = *myData[lvl];

      const DisjointBoxLayout& levelGrids = levelData.disjointBoxLayout();
      
      for (DataIterator dit = levelGrids.dataIterator(); dit.ok(); ++dit){
         const Box bx = levelGrids[dit()];
	 
         EBCellFAB& patchData       = levelData[dit()];
	 BaseFab<Real>& regularData = patchData.getSingleValuedFab(); // Regular data

	 // Iterate through all cells in bx
	 for (BoxIterator bit(bx); bit.ok(); ++bit){
	    const IntVect iv = bit();

	    regularData(iv, component) = ....
	 }

	 // Iterate through irregular cells later
	 for (VoFIterator vofit(...); vofit.ok(); ++vofit){
	    const VolIndex& vof = vofit();

	    patchData(vof, component)  = ...
	 }
      }
   }

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
   :width: 480px
   :align: center

   Example of stencils for computing gradients near embedded boundaries.
   The red stencil shows a regular 5-point stencil for computing the gradient on the coarse side of the refinement boundary; it reaches into the coarsened data beneath the fine level.
   The green stencil shows a similar 5-point stencil on the fine side of the refinement boundary; the stencil reaches over the refinement boundary and into one ghost cell.
   The blue stencils shows a much more complex stencil which is computed using a least squares reconstruction procedure. 

To compute gradients of a scalar, one can simply call ``AmrMesh::computeGradient(...)``.
See :ref:`Chap:AmrMesh` for details.   

Common data operations
----------------------

DataOps
_______

We have prototyped functions for many common data operations in a static class ``DataOps`` (see :file:`$DISCHARGE_HOME/Source/Utilities/CD_DataOps.H`).
For example, setting the value of a data holder can be done with

.. code-block:: c++

   EBAMRFluxData cellData;
   EBAMRFluxData fluxData;
   EBAMRIVData   irreData;
   
   DataOps::setValue(cellData, 0.0);
   DataOps::setValue(fluxData, 1.0);
   DataOps::setValue(irreData, 2.0);

Many functions are available in ``DataOps``.
Common data operations should always be put in this class. 

Copying data
____________

To copy data, one may use the ``EBAMRData<T>::copy(...)`` function *or* ``DataOps``.
These differ in the sense that ``DataOps`` will always do a local copy, and thus the data that is copied *must* be defined on the same realm.
Runtime errors will occur otherwise.
If you call ``EBAMRData<T>::copy(...)``, the data holders will first check if they are both defined on the same realm.
If they are, a purely local copy is perform.
Communication copies involving MPI are performed otherwise.

Coarsening
__________

Conservative coarsening data is done using the ``averageDown(...)`` function in ``AmrMesh``, with signatures as follows:

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

See :ref:`Chap:AmrMesh` for additional details. 

Ghost cells
___________

Filling ghost cells is done using the ``interpGhost(...)`` function in ``AmrMesh``.

.. code-block::

   // Interpolate ghost cells for multifluid cell-centered data.
   void interpGhost(MFAMRCellData& a_data, const std::string a_realm) const;

   // Interpolate ghost cells for multifluid cell-centered data   
   void interpGhost(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase) const;

This will fill the specified number of ghost cells using data from the coarse level only, using piecewise linear interpolation (with limiters).

As an alternative, one can interpolate a single layer of ghost cells using the multigrid interpolator (see :ref:`Chap:MultigridInterpolation`).
In this case only a single layer of ghost cells are filled in regular regions, but additional ghost cells (up to some specified range) are filled near the EB.
This is often required when computing gradients (to avoid reaching into invalid cut-cells).
See :ref:`Chap:Gradients`.
