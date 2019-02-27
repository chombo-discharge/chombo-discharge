.. _Chap:Tagging:

Cell tagging
------------

The :ref:`Chap:cell_tagger` class is too general for plug and play use since most users will not know how to interface into the solvers in order to fetch data. Because of this, we've added two simpler-to-use interfaces below this.

.. toctree::
   :maxdepth: 1

   UltralwTagger
   FullTagger

The ``ultralw_tagger`` is an abstract class where one can tag cells based on refinement criteria derived from the electric field only. The ``full_tagger`` class has access to much more solver data and can be used for more general cell tagging. Our favorite cell taggers are discussed below. 
   
streamer_tagger
_______________

The ``streamer_tagger`` is derived directly from ``cell_tagger`` and uses two different refinement criteria. One criteria is based on resolving gradients in the electric field by evaluation of

.. math::
   
   c_1 = \frac{\left|\nabla\left|\mathbf{E}\right|\right|}{\mathbf{E}}

and the other criteria uses
   
.. math::
      
   c_2 = \frac{S_e}{n_e\mathbf{v}_e}

If :math:`c_1\Delta x > \epsilon_1` or :math:`c_2\Delta x \epsilon_2`, the cell is refined. By default, non-refined cells are always coarsened. The user may select different :math:`\epsilon`'s for each level. 

field_tagger
____________

Currently, we mostly use field-based refinement and coarsening thresholds. The user may see the :doxy:`Doxygen API <classfield__tagger>` for details on this class. This class is a full implementation of :ref:`Chap:ultralw_tagger` and provides input parameters for controlling how refinement and coarsening is done.

.. literalinclude:: links/field_tagger.options

The *field_tagger* class uses curvature and relative-magnitude based refinement criteria. This class defines a single tracer field :math:`T = |\mathbf{E}|/\textrm{max}\left(|\mathbf{E}|\right)` and examines the curvature

.. math::
   
   c = \frac{\left|\nabla T\right|\Delta x}{T}

and the relative magnitude

.. math::

   r = T

If :math:`c` *or* :math:`r` are greater than ``refine_curvature`` or ``refine_magnitude``, a cell is refined. Likewise, if :math:`c` *and* :math:`r` are smaller than ``coarsen_curvature`` and ``coarsen_magnitude``, a cell is coarsened. 
