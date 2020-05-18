.. _Chap:MeshData:

Understanding mesh data
=======================

Most datastructures in `PlasmaC` are typedef'ed datastructures from Chombo, although a few new ones are also included.
For mesh data, the most central datastructure is the ``EBAMRCellFAB``, which holds cell-centered data across an AMR hierarchy which may or may not include embedded boundaries.

This class is typedef'ed as

.. code-block:: c++

   typedef Vector<RefCountedPtr<LevelData<EBCellFAB> > >

where the ``Vector`` holds a set of data on each AMR level; the data is allocated with a smart pointer called ``RefCountedPtr`` which points to a ``LevelData`` template structure, see :ref:`Chap:Basics`. 
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
