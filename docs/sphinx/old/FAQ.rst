.. _Chap:FAQ:

FAQ
===

This is a list to solutions of some common problems:

.. rubric:: How do I...

... change the spatial discretization?
   Please see the :ref:`Chap:amr_mesh` chapter.

... adjust the permittivity?
   The permittivity is a part of your :ref:`Chap:computational_geometry`. If your geometry supports dielectrics, your should be able to adjust it from there.

... restart a simulation?
   This is possible by writing checkpoint files and having :ref:`Chap:plasma_engine` start from one of your checkpoint filese. . Please see :ref:`Chap:RestartingSimulations` for further details. PS: There are *some* limitations to what parameters your change during your restarts.

... add another AMR level during restarts?
   You can't do this *unless* your finest level is lower than ``amr_mesh.max_amr_depth``. Please see :ref:`Chap:amr_mesh` for further details. 

... make my simulations run faster?
   Depending on your problem, you can try a few things. Try experimenting with different blocking factors by adjusting ``amr_mesh.blocking_factor`` and ``amr_mesh.max_box_size``. Please see :ref:`Chap:amr_mesh` to see what these parameters do. You can also try to manually remove some mesh by using :ref:`Chap:geo_coarsener`. If the convergence rate for the Poisson equation is low (:math:`<5`) then you can try experimenting with the level at which geometric multigrid enters the bottom solver. See :ref:`Chap:PoissonSolvers` for details.

... implement a new geometry?
   This is the work of :ref:`Chap:computational_geometry`. You can also check out the :ref:`Chap:NewSimulations` chapter to see how things are done.

... implement or alter a plasma kinetic scheme?
   This is the reponsibility of :ref:`Chap:plasma_kinetics`. You should also check out the worked examples in the :ref:`Chap:NewSimulations` chapter. 

... change the units in PlasmaC?
   You can't. We use SI units all over the place.

... speed up the geometry generation, or reduce the memory footprint?
   There are various things you can try. Try using a smaller ``amr_mesh.max_ebis_box_size`` or using the recursive box division algorithm for the index space generation. 

.. rubric:: Why does...
	      
... geometry generation take forever, or even crash?
   You probably have too many AMR levels for your core count. There are memory limitations to the maximum AMR depth that can be used for embedded boundary applications. See :ref:`Chap:EBMesh` for details. 

