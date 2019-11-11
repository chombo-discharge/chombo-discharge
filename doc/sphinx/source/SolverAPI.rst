.. _Chap:SolverAPI:

Solver API
==========

Here, we discuss the API of the base solvers on top of which `PlasmaC` is implemented. We only consider the public member functions that may appears through the `PlasmaC` core interface so that the core functions in `PlasmaC` can be understood for non-developers.


Minimal data structure guide
----------------------------

Most datastructures in `PlasmaC` are typedef'ed datastructures from Chombo. For example, the most central datastructure is the ``EBAMRCellFAB``, which holds cell-centered data across an AMR hierarchy which may or may not include embedded boundaries. This class is typedef'ed as

.. code-block:: c++

   typedef Vector<RefCountedPtr<LevelData<EBCellFAB> > >

where the ``Vector`` holds a set of data on each AMR level; the data is allocated with a smart pointer called ``RefCountedPtr`` which points to a ``LevelData`` template structure. The ``LevelData`` structure holds data on each grid patch on a specific AMR level. For example, in this case  ``LevelData<EBCellFAB>`` holds an ``EBCellFAB`` at each grid patch. The ``EBCellFAB`` is essentially a 2D or 3D array with cell-centered data (with appropriate index offsets). So, in order to retrieve data at a particular patch and grid level, one would perform the following:

.. code-block:: c++

   // Assume ebamrdata is an EBAMRCellFAB object
   EBCellFAB& data = (*ebamrdata[lvl])[dataIndex],

where ``dataIndex`` is an index used to retrieve data on a particular grid patch in ``LevelData``. Note that this will cause an error if a rank tries to retrieve data that it does not own. You will only see the above data retrieval inside iterators that iterate over the patches that are owned by specific MPI ranks. 

Other important datastructures are typedef'ed in the same way:

.. code-block:: c++

   typedef Vector<RefCountedPtr<LevelData<EBCellFAB> > >        EBAMRCellData; // Cell-centered EBAMR data
   typedef Vector<RefCountedPtr<LevelData<EBFluxFAB> > >        EBAMRFluxData; // Face-centered EBAMR data
   typedef Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > EBAMRIVData;   // Irregular cell EBAMR data
   typedef Vector<RefCountedPtr<LevelData<DomainFluxIFFAB > > > EBAMRIFData;   // Irregular face EBAMR data
   typedef Vector<RefCountedPtr<LevelData<MFCellFAB> > >        MFAMRCellData; // Like EBAMRCellData, but for two phases
   typedef Vector<RefCountedPtr<LevelData<MFFluxFAB> > >        MFAMRFluxData; // Like EBAMRFluxData, but for two phases
   typedef Vector<RefCountedPtr<LevelData<MFBaseIVFAB> > >      MFAMRIVData;   // Like EBAMRIVData, but for two phases

There are many more data structures in place, but the above data structures are the most commonly used ones. 
