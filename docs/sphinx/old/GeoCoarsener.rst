.. _Chap:geo_coarsener:

geo_coarsener
-------------

:ref:`Chap:geo_coarsener` is a plug-and-play class for removing "geometric" tags that were fetched by :ref:`Chap:plasma_engine`. This class is very simple: The user defines a number of boxes in space where boundary tags are removed down to a specified level. In this way, the user may remove tags in regions where the solution is well-behaved. The interface to :ref:`Chap:geo_coarsener` is

.. literalinclude:: links/geo_coarsener.options
		    
In the above, the user may modify ``num_boxes`` in order to specify a number of boxes where tags will be removed. As input, the boxes take the low and high corners, and the specified level on which tags are removed. Note that if tags are generated on level :math:`n`, the grid will have depth :math:`n+1`. Thus, if the user wishes to remove *all* boundary tags in a box, he must set ``boxN_lvl`` to zero. 
