.. _Chap:DesignOverview:

Overview
========

The major design components in ``chombo-discharge`` are:

#. :ref:`Chap:Driver` for running simulations.
#. :ref:`Chap:AmrMesh` for encapsulating (almost) all AMR and EB functionality in a common core class.
#. :ref:`Chap:TimeStepper` for integrating the equations of motion.
#. :ref:`Chap:ComputationalGeometry` for representing computational geometries (such as electrodes and dielectrics).
#. :ref:`Chap:CellTagger` for flagging cells for refinement and coarsening.
#. :ref:`Chap:GeoCoarsener` for manually removing refinement flags along the EB surface.

.. caution::

   In the current version of ``chombo-discharge``, :ref:`Chap:GeoCoarsener` has become redundant (due to improvements in the algorithms).
   It will be removed in future versions. 

