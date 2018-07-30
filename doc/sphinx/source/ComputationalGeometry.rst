.. _Chap:computational_geometry:

computational_geometry
----------------------

:ref:`Chap:computational_geometry` is the class that implements that geometry. In PlasmaC, we use level-set functions for description of surfaces. Please refer to :ref:`Chap:MiniApplications` for descriptions on how to implement new geometries.

:ref:`Chap:computational_geometry` is not an abstract class; if you pass in an instance of :ref:`Chap:computational_geometry` (rather than a casted instance), you will get a regular geometry without embedded boundaries.

By default, there are no input options available for :ref:`Chap:computational_geometry`, although inherited classes that actual implement a non-regular geometry will typically have many. Please see :ref:`Chap:MiniApplications` for further information. 
