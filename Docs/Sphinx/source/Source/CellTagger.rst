.. _Chap:CellTagger:

CellTagger
===========

The ``CellTagger`` class is responsible for flagging grid cells for refinement or coarsening.
If the user wants to implement a new refinement or coarsening routine, he will do so by writing a new derived class from ``CellTagger``.
The ``CellTagger`` parent class is a stand-alone class - it does not have a view of ``AmrMesh``, ``Driver``, or ``TimeStepper``.
Since refinement is intended to be quite general, the user is responsible for providing ``CellTagger`` with the appropriate depedencies.

.. tip::

   `CellTagger C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classCellTagger.html>`_

Refinement flags live in a data holder called ``EBAMRTags`` inside of ``Driver``.
This data is typedef'ed as

.. code-block:: c++

   typedef Vector<RefCountedPtr<LayoutData<DenseIntVectSet> > > EBAMRTags;

For performance reasons, ``DenseIntVectSet`` onlys store refinement cells on a per-patch basis.
It is not possible to add a grid cell to a ``DenseIntVectSet`` if it falls outside the grid patch.
Furthermore, the ``EBAMRTags`` structure is a distributed data structure which makes sure that each MPI rank is aware of the refinement flags on the corresponding grid patches.
The flags themselves are owned by :ref:`Chap:Driver`, and they are copied onto the new grids in the regrid step. 

User interface
--------------

To implement a new ``CellTagger``, the following functions must be implemented:

.. code-block:: c++

  virtual void regrid() = 0;
  virtual void parseOptions() = 0;
  virtual void parseRuntimeOptions();
  virtual bool tagCells(EBAMRTags& a_tags) = 0;

Users can also parse run-time options and even have ``CellTagger`` write to plot files, by implementing

.. code-block:: c++

  virtual int getNumberOfPlotVariables();
  virtual void writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp);

This is primarily useful for debugging the tracer fields that are used for flagging cells for refinement.

tagCells
________

When the regrid routine enters, the ``CellTagger`` will be asked to generate the refinement flags through a function

.. code-block:: c++

   bool tagCells(EBAMRTags& a_tags) = 0;

This routine should add cells that will be refined, and remove cells that will be coarsened.

regrid
______

``regrid`` is called by :ref:`Chap:Driver` during regrids.
The existence of this routine is due to the common usage pattern where ``CellTagger`` holds some auxiliary mesh data that is used when evaluating refinement and coarsening criteria.
This routine should reallocate such temporary storage during regrids.

parseOptions
____________

``parseOptions`` is called by :ref:`Chap:Driver` when setting the ``CellTagger``.
This routine should parse options into the ``CellTagger`` instance.

parseRuntimeOptions
___________________

``parseRuntimeOptions`` is called by :ref:`Chap:Driver` after each grid step.
This routine will re-read class options into the ``CellTagger`` instance.
See :ref:`Chap:RuntimeConfig` for further details.

getNumberOfPlotVariables
________________________

``getNumberOfPlotVariables`` will return the number of plot variables that ``CellTagger`` will write to plot files. 

writePlotData
_____________

``writePlotData`` will write the plot data to the provided data holder.
The functionality is the same as for :ref:`Chap:TimeStepper`. 

Restrict tagging
----------------

It is possible to prevent ``CellTagger`` from adding refinement flags in specified regions. 
The default behavior is to add a number of boxes where refinement and coarsening is allowed:

.. code-block:: bash

   MyCellTagger.num_boxes   = 0            # Number of allowed tag boxes (0 = tags allowe everywhere)
   MyCellTagger.box1_lo     = 0.0 0.0 0.0  # Only allow tags that fall between
   MyCellTagger.box1_hi     = 1.0 1.0 1.0  # these two corners

Here, ``MyCellTagger`` is a placeholder for the name of the class that is used.
By adding restrictive boxes, tagging will only be allowed inside the specified box corners ``box1_lo`` and ``box1_hi``.
More boxes can be specified by following the same convention, e.g. ``box2_lo`` and ``box2_hi`` etc.

Adding a buffer
---------------

By default, each MPI rank can only tag grid cells where it owns data.
This has been done for performance and communication reasons.
Under the hood, the ``DenseIntVectSet`` is an array of boolean values on a patch which is very fast and simple to communicate with MPI. 
Adding a grid cell for refinement which lies outside the patch will lead to memory corruptions.
It is nonetheless still possible to do this by growing the final generated tags like so:

.. code-block:: bash
		
   MyCellTagger.buffer = 4 # Add a buffer region around the tagged cells

Just before passing the flags into ``AmrMesh`` grid generation routines, the tagged cells are put in a different data holder (``IntVectSet``) and this data holder *can* contain cells that are outside the patch boundaries. 
