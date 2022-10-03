.. _Chap:GeoCoarsener:

GeoCoarsener
============

.. warning::

   ``GeoCoarsener`` is a relic of an early version of ``chombo-discharge``.
   It will be completely removed in a future version. 

The ``GeoCoarsener`` class permits manual removement of refinement flags along a geometric surface.
To remove these so-called *geometric tags*, the user specifies boxes in space where geometric tags will be removed:

.. code-block:: text
		
   GeoCoarsener.num_boxes   = 1            # Number of coarsening boxes (0 = don't coarsen)
   GeoCoarsener.box1_lo     = -1 -1 -1     # Lower-left corner 
   GeoCoarsener.box1_hi     =  1  1  1     # Upper-right corner
   GeoCoarsener.box1_lvl    = 0            # Remove tags down to this level. 
   GeoCoarsener.box1_inv    = false        # Flip removal.

If users want more boxes, they can specify it using the same syntax, e.g.

.. code-block:: text

   GeoCoarsener.num_boxes   = 2            # Number of coarsening boxes (0 = don't coarsen)
   
   GeoCoarsener.box1_lo     = -1 -1 -1     # Lower-left corner 
   GeoCoarsener.box1_hi     =  1  1  1     # Upper-right corner
   GeoCoarsener.box1_lvl    = 0            # Remove tags down to this level. 
   GeoCoarsener.box1_inv    = false        # Flip removal.

   GeoCoarsener.box2_lo     =  2  2  2     # Lower-left corner 
   GeoCoarsener.box2_hi     =  3  3  3     # Upper-right corner
   GeoCoarsener.box2_lvl    = 0            # Remove tags down to this level. 
   GeoCoarsener.box2_inv    = false        # Flip removal.		
