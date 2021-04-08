.. _Chap:SigmaSolver:

Surface charge solver
=====================

In order to conserve charge on solid insulators, `PlasmaC` has a solver that is defined on the gas-dielectric interface where the surface charge is updated with the incoming flux

.. math::
   F_\sigma(\phi) = \sum_{\phi}q_\phi F_{\textrm{EB}}(\phi),

where :math:`q_\phi` is the charge of a species :math:`\phi`. This ensures strong conservation on insulating surfaces.
