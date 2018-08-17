.. _Chap:Bugs:

Bugs
====

We don't know if PlasmaC is error-free (probably not). However, we take pains to remove all bugs that we find. If you encounter cases that don't work like you expect, we recommend that you :ref:`Chap:Contact`. 

Below, we present a list of some currently known bugs or missing features. 

... Sometimes when I place dielectrics, the Poisson solver doesn't produce the correct solution.
   In PlasmaC, dielectric interface cells are identified by the intersecting two of Chombo`s EBIndexSpace objects, which fails if the dielectric surface falls *exactly* on a grid face. We are working on an internal fix that handles this case. In the meantime, there is an easy workaround; simply displace the level set function away from the surface by a very small amount. 
