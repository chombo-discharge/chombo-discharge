.. _Chap:Bugs:

Bugs
====

We don't know if PlasmaC is error-free (probably not). However, we take pains to remove all bugs that we find. If you encounter cases that don't work like you expect, we recommend that you :ref:`Chap:Contact`. 

Below, we present a list of some currently known bugs, missing features, or performance tips. 

... Sometimes when I place dielectrics, the Poisson solver doesn't produce the correct solution.
   In PlasmaC, dielectric interface cells are identified by the intersecting two of Chombo`s EBIndexSpace objects, which fails if the dielectric surface falls *exactly* on a grid face. We are working on an internal fix that handles this case. In the meantime, there is an easy workaround; simply displace the level set function away from the surface by a very small amount.

... Sometimes, the Poisson solver doesn't converge.
   Geometric multigrid with embedded boundaries is a dark art. There are two things one should be aware of:

* In PlasmaC, the lower levels of multigrid are created by coarsening of the coarsest AMR level. How deep you can coarsen generally depends on your coarse-level decomposition (i.e. trough ``amr.max_box_size`` and ``amr.coarsest_domain``), and specification of multigrid coarsening through the specific solvers. In ``amr_mesh`` there is an option ``mg_coarsen`` that allows you to pre-coarsen some levels that will be used in multigrid. This works by coarsening the ``amr.coarsest_domain`` option a number of times by decomposing the coarsened domain with the ``amr.max_box_size`` setting. In principle, by setting this to a large value your lowest level multigrid domain can be equal to ``amr.max_box_size``. In many cases this is not enough, and you may want even more multigrid levels in order to achieve better convergence. This is done by adjusting the bottom solver depth in the elliptic solvers. This process is different from above since the lowest multigrid level may need to be less than ``amr.blocking_factor``. To do this, the lowest levels are coarsened by coarsening the grids directly, without changing the processor ownership. For example, if you have ``amr.coarsest_domain = 1024 1024`` and ``amr.mg_coarsen = 4``, the lower level grid coarsening will occur at a domain which is :math:`64\times 64`, which for ``amr.max_box_size = 32`` implies a processor grid which is :math:`2\times 2`. The lowest levels of multigrid will retain the processor grid, but the :math:`64\times64` domain is first coarsened to a :math:`16\times16` domain, then into a :math:`8\times8` domain, and so on. Often, performance of the Poisson solver is quite sensitve to these settings, and you should typically aim at getting your bottom solver down to e.g. an :math:`8\times 8` domain. 


  
* We do not support arbitrary grid configurations. For robust performance, you should refine all the embedded boundaries down to the finest level. This might not always be an option, and you might see degraded performance for the multigrid solvers when you have EBCF crossings. 



