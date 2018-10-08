.. _Chap:Bugs:

Bugs
====

We don't know if PlasmaC is error-free (probably not). However, we take pains to remove all bugs that we find. If you encounter cases that don't work like you expect, we recommend that you :ref:`Chap:Contact`. 

Below, we present a list of some currently known bugs, missing features, or performance tips. 

... Sometimes when I place dielectrics, the Poisson solver doesn't produce the correct solution.
   In PlasmaC, dielectric interface cells are identified by the intersecting two of Chombo`s EBIndexSpace objects, which fails if the dielectric surface falls *exactly* on a grid face. We are working on an internal fix that handles this case. In the meantime, there is an easy workaround; simply displace the level set function away from the surface by a very small amount.

... Sometimes, the Poisson solver doesn't converge.
   Geometric multigrid with embedded boundaries is a dark art. There are two things one should be aware of:

   * In PlasmaC, the lower levels of multigrid are created by coarsening of the coarsest AMR level. How deep you can coarsen generally depends on your coarse-level decomposition (i.e. trough ``amr.max_box_size`` and ``amr.coarsest_domain``). For example, a box size of 16 on the coarsest level can be coarsened 4 times to a box size of 1. This implies that a base domain of :math:`128^3` cells can be coarsened to an :math:`8^3` domain. Likewise, if you had a box size of :math:`32`, you could coarsen down to :math:`4^3`. Typically, the performance of geometric multigrid is better the further you can coarsen.
  
   * We do not support arbitrary grid configurations. For robust performance, you should refine all the embedded boundaries down to the finest level. This might not always be an option, and you might see degraded performance for the multigrid solvers when you have EBCF crossings. 
