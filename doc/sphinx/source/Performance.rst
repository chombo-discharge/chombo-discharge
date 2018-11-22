.. _Chap:Performance:

Performance
-----------

The performance of ``PlasmaC`` can be widely varying for different plasma kinetics and different geometries. Under the hood of ``PlasmaC``, there is much more going on than one might actual think. To demonstrate this in detail, the list below shows (almost) all the steps for a single Euler advance of the fluid equations, coupled with the Poisson and RTE equations.

For advancing the solver states from :math:`t^k` to :math:`t^{k+1}`, we do the following

1. Compute the electric field.
   
   a. Compute the cell-centered electric field :math:`\mathbf{E}` using centered differences :math:`E_i = -\frac{\phi_{i+1}-\phi_{i-1}}{2\Delta x}`.
   b. Compute the `face-centered` field by averaging the computed cell-centered field onto face centers. On domain boundaries where there are not enough cells available, we use second order backward/forward differences.
   c. Obtain the electric field at the embedded boundary centroids in order to handle charge injection. Do this by interpolating the cell-centered :math:`\mathbf{E}` to embedded boundary centroids. 

