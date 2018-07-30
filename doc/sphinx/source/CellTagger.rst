.. _Chap:cell_tagger:

cell_tagger
-----------

The :ref:`Chap:cell_tagger` class handles tagging of cells inside the gas phase of the simulation region. Currently, we do not support tagging of cells inside the solid phase, although this would be straightforward to implement. If this feature is desired, please contact us (see :ref:`Chap:contact` for contact information).

In PlasmaC, :ref:`Chap:cell_tagger` is an abstract class that the user must implement if he wishes to change how cells are tagged. The user must defined the number of desired tracer fields (through the constructor) and then implement three functions:


.. code-block:: c++

  /*!
    @brief Compute tracer fields
  */
  virtual void compute_tracers() = 0;
		
  /*!
    @brief Coarsen a cell based on a tracer field
  */
  virtual bool coarsen_cell(const RealVect&         a_pos,
			    const Real&             a_time,
			    const Real&             a_dx,
			    const int&              a_lvl,
			    const Vector<Real>&     a_tracer,
			    const Vector<RealVect>& a_grad_tracer) = 0;

  /*!
    @brief Refine a cell based on a tracer field
  */
  virtual bool refine_cell(const RealVect&         a_pos,
			   const Real&             a_time,
			   const Real&             a_dx,
			   const int&              a_lvl,
			   const Vector<Real>&     a_tracer,
			   const Vector<RealVect>& a_grad_tracer) = 0;

Internally, the regrid happens in the following way: In the regrid stage, :ref:`Chap:plasma_engine` will perform a call to :ref:`Chap:cell_tagger` to update a number of *tracer* fields: This is done in the *compute_tracers* function above. The tracer fields are scalar fields that the user defines: A tracer field may for example be :math:`T = |\mathbf{E}|`. Once those fields are updated, :ref:`Chap:cell_tagger` will look through all cells in the domain and perform refinement and coarsening based on the implementations of *coarsen_cell* and *refine_cell*.

For example, if the user has updated a tracer :math:`T = |\mathbf{E}|`, an example implemention of *refine_cell* is

.. code-block:: c++
		
		bool refine_cell(const RealVect&         a_pos,
            		         const Real&             a_time,
				 const Real&             a_dx,
				 const int&              a_lvl,
				 const Vector<Real>&     a_tracer,
				 const Vector<RealVect>& a_grad_tracer) {

		   return (a_grad_tracer[0]*a_dx)/a_tracer[0] > 0.5) ? true : false;
		}

In this implementation, the function will examine the local curvature :math:`\left|\nabla|\mathbf{E}|\right|\Delta x/|\mathbf{E}|` and refine if it is larger than :math:`0.5`. 


For full flexibility, :ref:`Chap:cell_tagger` has been granted solver access so that the user may fetch solution data directly from the solvers. For most users, this is uneccesarily complicated. We therefore have some mid-level implementation classes that simplify this process by fetching some frequently-used data to more well-defined interfaces (as opposed to *compute_tracers*). Please see :ref:`Chap:Solvers` for details.

There are options in the :ref:`Chap:cell_tagger` base class that permits the user to restrict tagging only to certain spatial regions by modifying the following parameters:

.. literalinclude:: links/cell_tagger.options

In the above, the user may define an arbitrary number of boxes in which tagging is *allowed*. If you do not specify a box, i.e. if **num_boxes** is zero, tagging is allowed everywhere. If specify one or more boxes, the **boxN_lo** and **boxN_hi** parameters indicate the valid tagging regions. Note that the boxes are not level-specific, since this is controlled through *coarsen_cell* and *refine_cell*, respectively. 
