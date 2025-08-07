.. _Chap:Basics:

``Chombo-3`` basics
===================

To fully understand this documentation the user should be familiar with ``Chombo``.
This documentation uses class names from ``Chombo`` and the most relevant ``Chombo`` data structures are summarized here.
What follows is a *very* brief introduction to these data structures, for in-depth explanations please see the `Chombo manual <https://github.com/applied-numerical-algorithms-group-lbnl/Chombo_3.2/tree/master/doc>`_.

Real
----

``Real`` is a typedef'ed structure for holding a single floating point number.
Compiling with double precision will typedef ``Real`` as ``double``, otherwise it is typedef'ed as ``float``.

RealVect
--------

``RealVect`` is a spatial vector.
It holds two ``Real`` components in 2D and three ``Real`` components in 3D.
The ``RealVect`` class has floating point arithmetic, e.g., addition, subtraction, multiplication etc.

Most of ``chombo-discharge`` is written in dimension-independent code, and for cases where ``RealVect`` is initialized with components the constructor uses ``Chombo`` macros for expanding the correct number of arguments.
For example

.. code-block:: c++

   RealVect v(D_DECL(vx, vy, vz));

will expand to ``RealVect v(vx,vy)`` in 2D and ``RealVect v(vx, vy, vz)`` in 3D.

IntVect
-------

``IntVect`` is an integer spatial vector, and is used for indexing data structures. 
It works in much the same way as ``RealVect``, except that the components are integers.

Box
---

The ``Box`` object describes a box in Cartesian space.
The boxes are indexed by the low and high corners, both of which are an ``IntVect``.
The ``Box`` may be cell-centered or face-centered.
To turn a cell-centered ``Box`` into a face-centered box one would do

.. code-block:: c++

   Box bx(IntVect::Zero, IntVect::Unit); // Default constructor give cell centered boxes
   bx.surroundingNodes():                // Now a cell-centered box

This will increase the box dimensions by one in each coordinate direction, and set the centering to face centered.

``Box`` is frequently used throughout ``chombo-discharge`` when iterating over grid cells.

EBCellFAB and FArrayBox
-----------------------

The ``EBCellFAB`` object is an array for holding cell-centered data in an embedded boundary context.
The ``EBCellFAB`` has two data structures: An ``FArrayBox`` which is a Cartesian array, and a additional data structure that holds data in multi-valued grid cells.
Doing arithmetic with ``EBCellFAB`` usually requires one to iterate over all the cell in the ``FArrayBox``, and then also to iterate over the *irregular cells* (i.e., cut-cells) later.
A ``VoFIterator`` is an object designed to iterate over the irregular cells, and must be defined by the set of cells that it will iterate over and the graph that describes the cell connectivity.

.. important::

   The ``FArrayBox`` stores the data in column major order.

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

``RefCountedPtr<T>`` is a pointer class in ``Chombo`` with reference counting. 
That is, when objects that hold a reference to some ``RefCountedPtr<T>`` object goes out of scope the reference counter is decremented.
If the reference counter reaches zero, the object that ``RefCountedPtr<T>`` points to it is deallocated.
Using ``RefCountedPtr<T>`` is much preferred over using a raw pointer ``T*`` to 1) avoid memory leaks and 2) compress code since no explicit deallocations need to be called. 

.. tip::
   
   In modern C++-speak, ``RefCountedPtr<T>`` can be thought of as a *very* simple version of ``std::shared_ptr<T>``. 

DisjointBoxLayout
-----------------

The ``DisjointBoxLayout`` class describes a grid on an AMR level where all the boxes are *disjoint*, i.e., boxes which do not overlap.
``DisjointBoxLayout`` is built upon a union of non-overlapping boxes having the same grid resolution and with unique rank-to-box ownership.
The constructor is

.. code-block:: c++

   Vector<Box> boxes(...);  // Vector of disjoint boxes
   Vector<int> ranks(...);  // Ownership of each box
   
   DisjointBoxLayout dbl(boxes, ranks);

In simple terms,  ``DisjointBoxLayout`` is the decomposed grid on each level in which MPI ranks have unique ownership of specific parts of the grid.

The ``DisjointBoxLayout`` is not a distributed data structure, and each MPI rank knows about all the boxes and the box ownership on the entire AMR level.
However, ranks will only allocate data on the part of the grid that they own. 
Data iterators also exist, and the most common is to use iterators that only iterate over its own part of the ``DisjointBoxLayout``.

.. code-block:: c++

   DisjointBoxLayout dbl;
   for (DataIterator dit(dbl); dit.ok(); ++dit){
      // Do something
   }

Each MPI rank will then iterate *only* over the part of the grid where it has ownership.

A related data iterators is ``LayoutIterator``, which will iterate over all boxes in the grid:

.. code-block:: c++

   for (LayoutIterator lit = dbl.layoutIterator(); dit.ok(); ++dit){
      // Do something
   }

This is typically used if one wants to do some global operations, e.g., count the number of cells in the grid.
However, trying to use ``LayoutIterator`` to retrieve data that was allocated by different MPI rank is an error. 
   

LevelData
---------

The ``LevelData<T>`` template structure holds data on all the grid patches of one AMR level.
The data is distributed with the domain decomposition specified by ``DisjointBoxLayout``, and each patch contains exactly one instance of ``T``.
``LevelData<T>`` uses a factory pattern for creating the ``T`` objects, so if you have new data structures that should fit within ``LevelData<T>`` structure you must also implement a factory method for ``T``, as well as an appropriate linearization function for ``T``.

The ``LevelData<T>`` object provides the domain decomposition method in ``Chombo`` and ``chombo-discharge``.
Often, ``T`` is an ``EBCellFAB``.

To iterate over ``LevelData<T>`` one will use the data iterator above: 

.. code-block:: c++

   LevelData<T> myData;
   for (DataIterator dit(dbl); dit.ok(); ++dit){
      T& = myData[dit()];
   }

``LevelData<T>`` also includes the concept of ghost cells and exchange operations.


EBISLayout and EBISBox
----------------------

The ``EBISLayout`` holds the geometric information over one ``DisjointBoxLayout`` level.
Typically, the ``EBISLayout`` is used for fetching the geometric moments that are required for performing computations near cut-cells. 
``EBISLayout`` can be thought of as an object which provides all EB-related information on a specific grid level.
The EB information consists of, e.g., cell flags (i.e., is the cell a cut-cell?), volume fractions, normal vectors, etc. 
This information is stored in a class ``EBISBox``, which holds all the EB information for one specific grid patch.
To obtain the EB-information for a specific grid patch, one will call:

.. code-block:: c++

   EBISLayout ebisl;
   for (DataIterator dit(dbl); dit.ok(); ++dit){
      EBISBox& ebisbox = ebisl[dit()];
   }

where ``EBISBox`` contains the geometric information over only one grid patch.
One can thus think of the ``EBISLayout`` as a ``LayoutData<EBISBox>`` structure. 

As an example, to iterate over all the cut-cells defined for a cell-centered data holder an AMR-level one would do:

.. code-block:: c++

   constexpr int comp = 0;

   // Assume that these exist. 
   LevelData<EBCellFAB> myData;
   EBISLayout ebisl;

   // Iterate over all the patches on a grid level.
   for (DataIterator dit(dbl); dit.ok(); ++dit){
      const Box      cellBox = dbl[dit()];
      const EBISBox& ebisbox = ebisl [dit()];

      EBCellFAB& patchData = myData[dit()];      

      // Get all the cut-cells in the grid patch
      const IntVectSet& ivs = ebisbox.getIrregIVS(cellBox);
      const EBGraph&        = ebisbox.getEBGraph();

      // Define a VoFIterator for the cut-cells and iterate over all the cut-cells.
      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
         const VolIndex& vof = vofit();

	 patchData(vof, comp) = ...
      }
   }

Here, ``EBGraph`` is the graph that describes the connectivity of the cut cells.

.. _Chap:BaseIF:

BaseIF
------

The ``BaseIF`` is a ``Chombo`` class which encapsulates an implicit function (recall that all SDFs are also implicit functions, see :ref:`Chap:GeometryRepresentation`).
``BaseIF`` is therefore used for fundamentally constructing a geometric object.
Many examples of ``BaseIF`` are found in ``Chombo`` itself, and ``chombo-discharge`` includes additional ones.

To implement a new implicit function, the user must inherit from ``BaseIF`` and implement the pure function

.. code-block:: c++

   virtual Real BaseIF::value(const RealVect& a_point) const = 0;

The implemention should return a positive value if the point ``a_point`` is inside the object and a negative value otherwise.
