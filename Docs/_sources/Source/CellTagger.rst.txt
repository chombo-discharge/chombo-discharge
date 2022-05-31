.. _Chap:CellTagger:

CellTagger
===========

The ``CellTagger`` class is responsible for flagging grid cells for refinement or coarsening.
If the user wants to implement a new refinement routine, he will do so by writing a new derived class from ``CellTagger``.
The ``CellTagger`` parent class is a stand-alone class - it does not have a view of ``AmrMesh``, ``Driver``, or ``TimeStepper``.
Since refinement is intended to be quite general, the user is responsible for providing ``CellTagger`` with the appropriate depedencies.
For example, for streamer simulations we often use the electric field when tagging grid cells, in which case the user should supply either a reference to the electric field, the Poisson solver, or the time stepper.

Design
------
In ``chombo-discharge``, refinement flags live in a data holder called ``EBAMRTags`` inside of ``Driver``.
This data is typedef'ed as

.. code-block:: c++

   typedef Vector<RefCountedPtr<LayoutData<DenseIntVectSet> > > EBAMRTags;

The ``LayoutData<T>`` structure can be thought of as a ``LevelData<T>`` without possibilities for communication.
``CellTagger`` is an abstract class that the user *must* overwrite. 

When the regrid routine enters, the ``CellTagger`` will be asked to generate the refinement flags through a function

.. code-block:: c++

   bool tagCells(EBAMRTags& a_tags) = 0;

Typically, code for refinement looks something like:

.. code-block:: c++

   bool myCelltagger::tagCells(EBAMRTags& a_tags){

      for (int l 0; l <= finestLevel; lvl++){
         const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
	 
         for (DataIterator dit(dbl); dit.ok(); ++dit){
	     DenseIntVectSet& boxTags = (*a_tags[lvl])[dit()];
	     const Box& box           = dbl[dit()];

	     for (BoxIterator bit(box); bit.ok(); ++bit){
	        const IntVect iv = bit();

		const bool refineThisCell = myRefinementFunction(iv, ...);

		if(refineThisCell){
		   boxTags(iv) = true;
		}
		else{
		   boxTags(iv) = false;
		}
	     }
	 }

      return true;
   }

User interface
--------------

To implement a new ``CellTagger``, the following functions must be implemented:

.. code-block:: c++

  virtual void regrid() = 0;
  virtual void parseOptions() = 0;
  virtual void parseRuntimeOptions();
  virtual bool tagCells(EBAMRTags& a_tags) = 0;

Users can also parse run-time options and have ``CellTagger`` write to file, by implementing

.. code-block:: c++

  virtual int getNumberOfPlotVariables();
  virtual void writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp);

This is primarily useful for debugging the tracer fields that are used for flagging cells for refinement. 
   

Restrict tagging
----------------

It is possible to prevent cell tags in certain regions.
The default is simply to add a number of boxes where refinement and coarsening is allowed by specifying a number of boxes in the options file for the cell tagger:

.. code-block:: bash

   MyCellTagger.num_boxes   = 0            # Number of allowed tag boxes (0 = tags allowe everywhere)
   MyCellTagger.box1_lo     = 0.0 0.0 0.0  # Only allow tags that fall between
   MyCellTagger.box1_hi     = 0.0 0.0 0.0  # these two corners

Here, *MyCellTagger* is a placeholder for the name of the class that is used.
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
		
   MyCellTagger.buffer = 4 # Add a buffer region around the tagged cells

Just before passing the flags into ``AmrMesh`` grid generation routines, the tagged cells are put in a different data holder (``IntVectSet``) and this data holder *can* contain cells that are outside the patch boundaries. 
