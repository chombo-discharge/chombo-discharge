.. _Chap:MeshData:

Understanding mesh data
=======================

Most datastructures in `PlasmaC` are typedef'ed datastructures from Chombo, although a few new ones are also included.
For mesh data, the most central datastructure is the ``EBAMRCellFAB``, which holds cell-centered data across an AMR hierarchy which may or may not include embedded boundaries.

This class is typedef'ed as

.. code-block:: c++

   typedef Vector<RefCountedPtr<LevelData<EBCellFAB> > >

where the ``Vector`` holds a set of data on each AMR level; the data is allocated with a smart pointer called ``RefCountedPtr`` which points to a ``LevelData`` template structure.
In modern C++ speak, one may think of ``RefCountedPtr`` as a ``std::shared_ptr``.

The ``LevelData<T>`` template structure holds data on all the grid patches of one AMR level.
The data is distributed in parallel, and each patch contains exactly one instance of ``T``.
``LevelData`` is built upon a class called ``DisjointBoxLayout``, which is the collection of all grid patches one AMR level.
One may think of the ``DisjointBoxLayout`` as the grid on each level, and where MPI ranks take unique ownership of data on each patch. 

In the above typedef'ed example, ``LevelData<EBCellFAB>`` holds an ``EBCellFAB`` inside each grid patch. 
The ``EBCellFAB`` is essentially a 2D or 3D array with cell-centered data (with appropriate index offsets), and includes special handling in the case that a cut-cell is multivalued. 
So, in order to retrieve data at a particular patch and grid level, one would perform the following:

.. code-block:: c++

   // Assume ebamrdata is an EBAMRCellFAB object
   EBCellFAB& data = (*ebamrdata[lvl])[dataIndex],

where ``dataIndex`` is an index used to retrieve data on a particular grid patch in ``LevelData``.
Note that this will cause an error if a rank tries to retrieve data that it does not own.
You will only see the above data retrieval inside iterators that iterate over the patches that are owned by specific MPI ranks, like this:

.. code-block:: c++

   // Assume we have a LevelData<EBCellFAB> levelData built upon a DisjointBoxLayout called dbl
   for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      EBCellFAB& data = levelData[dit()];
   }

Here, the class ``DataIterator`` is a parallel iterator that iterates *only* over the patches owned by each MPI rank.

Other important datastructures are typedef'ed in much the same way as the cell-centered data:

.. code-block:: c++

   typedef Vector<RefCountedPtr<LevelData<EBCellFAB> > >        EBAMRCellData; // Cell-centered EBAMR data
   typedef Vector<RefCountedPtr<LevelData<EBFluxFAB> > >        EBAMRFluxData; // Face-centered EBAMR data
   typedef Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > EBAMRIVData;   // Irregular cell EBAMR data
   typedef Vector<RefCountedPtr<LevelData<DomainFluxIFFAB > > > EBAMRIFData;   // Irregular face EBAMR data
   typedef Vector<RefCountedPtr<LevelData<MFCellFAB> > >        MFAMRCellData; // Like EBAMRCellData, but for two phases
   typedef Vector<RefCountedPtr<LevelData<MFFluxFAB> > >        MFAMRFluxData; // Like EBAMRFluxData, but for two phases
   typedef Vector<RefCountedPtr<LevelData<MFBaseIVFAB> > >      MFAMRIVData;   // Like EBAMRIVData, but for two phases

There are many more data structures in place, but the above data structures are the most commonly used ones.
Here, ``EBAMRFluxData`` is precisely like ``EBAMRCellData``, except that the data is stored on *cell faces* rather than cell centers.
Likewise, ``EBAMRIVData`` is a typedef'ed data holder that holds data on each cut-cell center across the entire AMR hierachy.
In the same way, ``EBAMRIFData`` holds data on each face of all cut cells. 
