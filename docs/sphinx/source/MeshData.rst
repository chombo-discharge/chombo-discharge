.. _Chap:MeshData:

Understanding mesh data
=======================

Mesh datastructures in ``PlasmaC`` are derived from a class ``EBAMRData<T>`` which holds a typename ``T`` at every box at every level in the AMR hierarchy.
Internally, the data is stored as a ``Vector<RefCountedPtr<LevelData<T> > >`` and the user may fetch the data through ``EBAMRData<T>::get_data()``.
Here, the ``Vector`` holds a set of data on each AMR level; the data is allocated with a smart pointer called ``RefCountedPtr`` which points to a ``LevelData`` template structure, see :ref:`Chap:Basics`. 

The reason for having class encapsulation of mesh data is due to :ref:`Chap:realm`, so that we can only keep track on which realm the mesh data is defined.
Users will not have to interact with ``EBAMRData<T>`` directly, but will do so primarily through ``amr_mesh``.
``amr_mesh`` has functionality for defining most ``EBAMRData<T>`` types on one of the realms, and ``EBAMRData<T>`` itself it typically not used anywhere elsewhere within ``PlasmaC``.

A number of explicit template specifications also exist, and these are outlined below: 

.. code-block:: c++

   typedef EBAMRData<MFCellFAB>        MFAMRCellData;  // Cell-centered multifluid data
   typedef EBAMRData<MFFluxFAB>        MFAMRFluxData;  // Face-centered multifluid data
   typedef EBAMRData<MFBaseIVFAB>      MFAMRIVData;    // Irregular face multifluid data
   typedef EBAMRData<EBCellFAB>        EBAMRCellData;  // Cell-centered single-phase data
   typedef EBAMRData<EBFluxFAB>        EBAMRFluxData;  // Face-centered data in all coordinate direction
   typedef EBAMRData<EBFaceFAB>        EBAMRFaceData;  // Face-centered in a single coordinate direction
   typedef EBAMRData<BaseIVFAB<Real> > EBAMRIVData;    // Data on irregular data centroids
   typedef EBAMRData<DomainFluxIFFAB>  EBAMRIFData;    // Data on domain phases
   typedef EBAMRData<BaseFab<bool> >   EBAMRBool;      // For holding bool at every cell

There are many more data structures in place, but the above data structures are the most commonly used ones.
Here, ``EBAMRFluxData`` is precisely like ``EBAMRCellData``, except that the data is stored on *cell faces* rather than cell centers.
Likewise, ``EBAMRIVData`` is a typedef'ed data holder that holds data on each cut-cell center across the entire AMR hierachy.
In the same way, ``EBAMRIFData`` holds data on each face of all cut cells. 
