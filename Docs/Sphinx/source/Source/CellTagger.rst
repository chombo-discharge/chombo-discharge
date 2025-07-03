.. _Chap:CellTagger:

CellTagger
===========

The ``CellTagger`` class is responsible for flagging grid cells for refinement or coarsening, and is thus the main class that determines what the grids look like.
By default, ``CellTagger`` does not actually flag anything for refinement, so if the user wants to implement a new refinement or coarsening routine, one must do so by writing a new derived class from ``CellTagger``.
The ``CellTagger`` parent class is a stand-alone class - it does not have a view of :ref:`Chap:AmrMesh`, :ref:`Chap:Driver`, or :ref:`Chap:TimeStepper` and the user is responsible for providing ``CellTagger`` with these dependencies.

.. tip::

   Here is the `CellTagger C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classCellTagger.html>`_

Refinement flags live in a data holder called ``EBAMRTags``, which is a typedef:

.. code-block:: c++

   typedef Vector<RefCountedPtr<LayoutData<DenseIntVectSet> > > EBAMRTags;

See :ref:`Chap:Basics` for an explanation of how the individual templates.
Briefly, the outer vector in ``EBAMRTags`` indicates the grid level, whereas ``LayoutData`` is a data holder on each grid level and indexes the grid patches.
That is, ``LayoutData<DenseIntVectSet`` is a distribution of ``DenseIntVectSet`` on a level, whereas ``DenseIntVectSet`` is a data holder that stores the refinement flags on a grid patch.
For performance reasons, ``DenseIntVectSet`` onlys store refinement cells on a per-patch basis.
It is not possible to add a grid cell to a ``DenseIntVectSet`` if it falls outside the grid patch.
Note that ``EBAMRTags`` is *not* set up for communication, and so it is not possible to fill ghost cells regions from other patches.
The refinement flags themselves are owned by :ref:`Chap:Driver`, and are used to generate new grids that cover all the cell tags.

A part of the current C++ header file for ``CellTagger`` is included below, where we highlight the functions that must be implemented in order to create a new refinement method.

.. _CellTaggerPureFunctions:
.. literalinclude:: ../../../../Source/CellTagger/CD_CellTagger.H
   :language: c++
   :caption: Header file for ``CellTagger``.
   :lines: 49-62, 71-79
   :emphasize-lines: 6-7, 13-14, 22-23

Implementing a new ``CellTagger``
---------------------------------

To implement a new ``CellTagger``, the pure functions listed in :numref:`CellTaggerPureFunctions` must be implemented.
Below, we discuss these in turn:

regrid
______

``regrid`` is called by :ref:`Chap:Driver` during regrids.
The existence of this routine is due to the common usage pattern where ``CellTagger`` holds some auxiliary mesh (or particle) data that is used when evaluating refinement and coarsening criteria.
This routine should reallocate this  storage during regrids.

tagCells
________

When the regrid routine enters, the ``CellTagger`` will be asked to generate the refinement flags through the ``tagCells`` function in :numref:`CellTaggerPureFunctions`. 
This routine should add cells that will be refined, and remove cells that will be coarsened.
The return value of this function is a boolean that should return ``false`` if no new tags were found.
While it is possible to skip this test, we note that returning ``true`` also when no new refinement regions were found will trigger a full regrid.

parseOptions
____________

``parseOptions`` is called by :ref:`Chap:Driver` when setting the ``CellTagger``.
This routine should parse options into the ``CellTagger`` subclass.
Note that it is also useful to overwrite the function ``parseRuntimeOptions`` if one needs run-time configuration of the refinement criteria.

Plot data
_________

It is also possible to have ``CellTagger`` write data to plot files, and the interface for this is identical to the plot file interface in :ref:`Chap:TimeStepper`.
The three functions below must then be implemented:

.. _CellTaggerPureFunctions:
.. literalinclude:: ../../../../Source/CellTagger/CD_CellTagger.H
   :language: c++
   :lines: 87-108
   :emphasize-lines: 5-6, 11-12, 21-22

The interpretation of these functions is exactly the same as for :ref:`Chap:TimeStepper`, and we refer to the :ref:`Chap:TimeStepper` documentation for a detailed explanation.

Manual refinement and restriction
---------------------------------

The user can add manual refinement by specifying Cartesian spatial regions to be refined down to some grid level, by specifying the physical corners and the refinement level.
The input parameters in this case are

* ``num_ref_boxes`` for specifying how many such boxes will be parsed.
* ``ref_box<num>_lo`` and ``ref_box<num>_hi`` that determine the Cartesian region to be refined.
  Here, ``<num>`` is a placeholder for an integer.
  If the user specifies ``num_ref_boxes = 2``, then ``ref_box1_lo``, ``ref_box1_hi``, ``ref_box2_lo``, and ``ref_box2_hi`` must all be defined.
* ``ref_box<num>_lvl``, which specifies the refinement depth corresponding to box ``<num>``.

An example of the manual refinement syntax is given in :numref:`CellTaggerOptions`.

.. _CellTaggerOptions:
.. literalinclude:: ../../../../Source/CellTagger/CD_CellTagger.options
   :caption: Default class options for ``CellTagger``, including the manual refinement syntax and buffer region definition.
   :language: txt

It is possible to prevent ``CellTagger`` from adding refinement flags in specified regions by specifying ``num_tag_boxes``.
The default behavior is to add a number of boxes where refinement and coarsening is allowed, and prevent tags from being generated outside of these regions.
If no boxes are defined, tagging is allowed everywhere.
The syntax is the same as for ``num_ref_boxes``, e.g., one must define ``tag_box<num>_lo`` and ``tag_box<num>_hi``.
By adding restrictive boxes, tagging will only be allowed inside the specified box corners ``tag_box1_lo`` and ``tag_box1_hi``.
More boxes can be specified by following the same convention, e.g., ``tag_box2_lo`` and ``tag_box2_hi``.

Adding a buffer
---------------

By default, each MPI rank can only tag grid cells where it owns data, which has has been enforced due to performance and communication reasons.
Under the hood, the ``DenseIntVectSet`` is an array of boolean values on a patch which is very fast and simple to communicate with MPI. 
Adding a grid cell for refinement which lies outside the patch will lead to memory corruptions.
Frequently, however, it is necessary to add *buffer regions* to ensure that an area-of-interest around the tagged region is also refined.
This is done by growing the final generated tags by specifying the ``buffer``.
The syntax for this is given in :numref:`CellTaggerOptions`.
Just before passing the flags into ``AmrMesh`` grid generation routines, the tagged cells are put into a different data holder (``IntVectSet``), and this data holder *can* contain cells that are outside the patch boundaries.		     








