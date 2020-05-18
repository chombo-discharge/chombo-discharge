.. _Chap:cell_tagger:

cell_tagger
===========

The ``cell_tagger`` class is responsible for flagging grid cells for refinement or coarsening.
In ``PlasmaC``, refinement flags live in a data holder called ``EBAMRTags`` inside of ``driver``.
This data is typedef'ed as

.. code-block:: c++

   typedef Vector<RefCountedPtr<LayoutData<DenseIntVectSet> > > EBAMRTags;

The ``LayoutData<T>`` structure can be thought of as a ``LevelData<T>`` without possibilities for communication.
``cell_tagger`` is an abstract class that the user *must* overwrite - it is not possible to 

When the regrid routine enters, the ``cell_tagger`` will be asked to generate the refinement flags through a function

.. code-block:: c++

   bool tag_cells(EBAMRTags& a_tags) = 0;

If the user wants to implement a new refinement routine, he will do so by writing a new derived class from ``cell_tagger``.
The ``cell_tagger`` parent class is a stand-alone class - it does not have a view of ``amr_mesh``, ``driver``, or ``time_stepper``.
Since refinement is intended to be quite general, the user is responsible for providing ``cell_tagger`` with the appropriate depedencies.
For example, for streamer simulations we often use the electric field when tagging grid cells, in which case the user should supply either a reference to the electric field, the potential, the Poisson solver, or the plasma physics object. 

Restrict tagging
----------------

It is possible to prevent cell tags in certain regions.
The default is simply to add a number of boxes where refinement and coarsening is allowed by specifying a number of boxes in the options file for the cell tagger:

.. code-block:: bash

   my_cell_tagger.num_boxes   = 0            # Number of allowed tag boxes (0 = tags allowe everywhere)
   my_cell_tagger.box1_lo     = 0.0 0.0 0.0  # Only allow tags that fall between
   my_cell_tagger.box1_hi     = 0.0 0.0 0.0  # these two corners

Here, *my_cell_tagger* is a placeholder for the name of the class that is used.
By adding restrictive boxes, tagging will only be allowed inside the specified box corners ``box1_lo`` and ``box1_hi``.
More boxes can be specified by following the same convention, e.g. ``box2_lo`` and ``box2_hi`` etc.

Adding a buffer
---------------

By default, each MPI rank can only tag grid cells where he owns data.
This has been done for performance and communication reasons.
Under the hood, the ``DenseIntVectSet`` is an array of boolean values on a patch which is very fast and simple to communicate with MPI. 
Adding a grid cell for refinement which lies outside the patch will lead to memory corruptions.
It is nonetheless still possible to do this by growing the final generated tags like so:

.. code-block:: bash
		
   my_cell_tagger.buffer = 4 # Add a buffer region around the tagged cells

Just before passing the flags into ``amr_mesh`` grid generation routines, the tagged cells are put in a different data holder (``IntVectSet``) and this data holder *can* contain cells that are outside the patch boundaries. 
