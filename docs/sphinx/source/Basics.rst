.. _Chap:Basics:

Chombo basics
=============

To fully understand this documentation the user must be familiar with Chombo.
This documentation uses class names from Chombo and the most relevant Chombo data structures are summarized here.
What follows is a *very* brief introduction to these data structures, for in-depth explanations please see the Chombo manual. 

Real
----

``Real`` is a typedef'ed structure for holding a single floating point number.
Compilating with double precision will typedef ``Real`` as ``double``, otherwise it is typedef'ed as ``float`` 

RealVect
--------

``RealVect`` is a spatial vector.
It holds two ``Real`` components in 2D and three ``Real`` components in 3D.
The ``RealVect`` class has floating point arithmetic, e.g. addition, subtraction, multiplication etc.

Most of ``PlasmaC`` is written in dimension-independent code, and for cases where ``RealVect`` is initialized with components the constructor uses Chombo macros for expanding the correct number of arguments.
For example

.. code-block:: c++

   RealVect v(D_DECL(vx, vy, vz));

will expand to ``RealVect(vx,vy)`` in 2D and ``RealVect(vx, vy, vz)`` in 3D.

IntVect
-------

``IntVect`` is an integer spatial vector, and is used for indexing data structures. 
It works in much the same way as ``RealVect``, except that the components are integers.

Box
---

The ``Box`` object describes a Box in Cartesian space.
The boxes are indexed by the low and high corners, both of which are ``IntVect``s. 
The ``Box`` may be cell-centered or face-centered.
To turn a cell-centered ``Box`` into a face-centered box one would do

.. code-block:: c++

   Box bx(IntVect::Zero, IntVect::Unit); // Default constructor give cell centered boxes
   bx.surroundingNodes():                // Now a cell-centered box

This will increase the box dimensions by one in each coordinate direction.

EBCellFAB and FArrayBox
-----------------------

The ``EBCellFAB`` object is an array for holding cell-centered data in an embedded boundary context.
The ``EBCellFAB`` is defined over a ``Box`` and the cut-cell geometry in the same ``Box``. 
For cut-cells, there is a chance that the boundaries intersect the grid in such a way that a cell has more than one degree of freedom.
Data therefore needs to live on more complex data structures than simple arrays.

The ``EBCellFAB`` has two data structures: An ``FArrayBox`` that holds the data on the cell centers, and is supported by a graph that additionally holds data in cells that are multiply cut.
The ``FArrayBox`` is essentially a Fortran array that can be passed to Fortran for performance reasons.
Doing arithmetic with ``EBCellFAB`` usually requires one to iterate over all the cell in the ``FArrayBox``, and then to iterate over the *irregular cells* (i.e. cut-cells) later.
The ``VofIterator`` can iterate over any number cells, but is typically used only to iterate over the cut-cells in a ``Box``.
Typically, code for doing anything with the ``EBCellFAB`` looks like this:

.. code-block:: c++

   // Call fortran code
   FORT_DO_SOMETHING(....)

   // Iterate over cut-cells
   for (VoFIterator vofit(...); vofit.ok(); ++vofit){
      (...)
   }

Vector
------

``Vector<T>`` is a one-dimensional array with constant-time random access and range checking.
It uses ``std::vector`` under the hood and can access the most commonly used ``std::vector`` functionality through the public member functions.
E.g. to obtain an element in the vector

.. code-block:: c++

   Vector<T> my_vector(10, T());

   T& element = my_vector[5];

Likewise, ``push_back``, ``resize`` etc works in much the same way as for ``std::vector``.

RefCountedPtr
-------------

``RefCountedPtr<T>`` is a pointer class in Chombo with reference counting.
That is, when objects that hold a reference to some ``RefCountedPtr<T>`` object goes out of scope the reference counter is decremented.
If the reference counter reaches zero, the object that ``RefCountedPtr<T>`` points to it deallocated.
Using ``RefCountedPtr<T>`` is much preferred over using a bare pointer ``T*`` to 1) avoid memory leaks and 2) compress code since no explicit deallocations need to be called. 

In modern C++ speak, ``RefCountedPtr<T>`` can be thought of as a bare-bones ``std::shared_ptr`` (without the move semantics and so on). 

DisjointBoxLayout
-----------------

The ``DisjointBoxLayout`` describes a grid on an AMR level where all the boxes are *disjoint*, i.e. they don't overlap.
``DisjointBoxLayout`` is built upon a collection of boxes and the MPI rank ownership of those boxes.
The constructor is

.. code-block:: c++

   Vector<Box> boxes(...);  // Vector of disjoint boxes
   Vector<int> ranks(...);  // Ownership of each box
   
   DisjointBoxLayout dbl(boxes, ranks);

In simple terms,  ``DisjointBoxLayout`` is the decomposed grid on each level in which MPI ranks have unique ownership of specific parts of the grid.

The ``DisjointBoxLayout`` view is global, i.e. each MPI rank knows about all the boxes and the box ownership on the entire AMR level.
However, ranks will only allocate data on the part of the grid that they own. 
Data iterators also exist, and the most common is to use iterators that only iterate over the part of the ``DisjointBoxLayout`` that the specific MPI ranks own:

.. code-block:: c++

   for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      // Do something
   }

Each MPI rank will then iterate *only* over the part of the grid where it has ownership.

Other data iterators exist that iterate over the *whole* grid:

.. code-block:: c++

   for (LayoutIterator dit = dbl.layoutIterator(); dit.ok(); ++dit){
      // Do something
   }

This is typically used if one wants to do some global operation, e.g. count the number of cells in the grid or somesuch.
If you try to use ``LayoutIterator`` to retrieve data that was allocated locally, you will get a memory corruption. 
   

LevelData
---------

The ``LevelData<T>`` template structure holds data on all the grid patches of one AMR level.
The data is distributed with the domain decomposition specified by ``DisjointBoxLayout``, and each patch contains exactly one instance of ``T``.
``LevelData<T>`` uses a factory pattern for creating the ``T`` objects, so if you have new data structures that should fit the in ``LevelData<T>`` structure you must also implement a factory method for ``T``. 

To iterate over ``LevelData<T>`` one will use the data iterator above: 

.. code-block:: c++

   LevelData<T> myData;
   for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const DataIndex& d = dit();
      T& = myData[d];
   }

where ``[DataIndex]`` is an indexing operator for ``LevelData``.

``LevelData<T>`` also includes the concept of ghost cells and exchange operations.
Specifying ghost cells is primarily controlled in input scripts for simulations.




EBISLayout and EBISBox
----------------------

The ``EBISLayout`` holds the geometric information over one ``DisjointBoxLayout`` level.
Typically, the ``EBISLayout`` is used to obtain information about the cut-cells.
``EBISLayout`` also has an indexing operator that can be used to extract a pointer to the geometric information in only one of the boxes on the level.
This indexing operator is

.. code-block:: c++

   EBISLayout ebisl(...);
   for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const DataIndex& d = dit();

      EBISBox& ebisbox = ebisl[d];
   }

where ``EBISBox`` contains the geometric information over only one ``Box``. 

As an example, to iterate over all the cut-cells defined for a cell-centered data holder an AMR-level one would do:

.. code-block:: c++

   const int comp = 0;
   
   LevelData<EBCellFAB> myData(...);
   EBISLayout ebisl(...);
   
   for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const DataIndex& d = dit();
      const Box bx       = dbl.get(d);

      EBCellFAB& ebcell = myData[d];
      EBISBox& ebisbox  = ebisl[d];

      const IntVectSet& ivs = ebisbox.getIrregIVS(box);
      const EBGraph&        = ebisbox.getEBGraph();
      
      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
         const VolIndex& vof = vofit();

	 ebcell(vof, comp) = ...
      }
   }

Here, ``EBGraph`` is the graph that describes the connectivity of the cut cells. 
