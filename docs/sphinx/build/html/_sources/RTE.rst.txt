.. _Chap:RadiativeTransfer:

Radiative transfer
==================

Radiative transfer is supported in the diffusion (i.e. Eddington or Helmholtz) approximation and with Monte Carlo sampling of discrete photons.
The solvers share a common interface but since diffusion RTE is deterministic and discrete Monte Carlo photons are stochastic, not all temporal integration methods will support both.
The diffusion approximation relies on solving an elliptic equation in the stationary case and a parabolic equation in the time-dependent case, while the Monte-Carlo approach solves solves for fully transient or ''stationary'' transport. 

Diffusion approximation
-----------------------

In the diffusion approximation, the radiative transport equation is

.. math::

      \partial_t\Psi + \kappa\Psi - \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi\right) = \frac{\eta}{c},

which is called the Eddington approximation. The radiative flux is :math:`F = -\frac{c}{3\kappa}\nabla \Psi`.
We do not currently support flux-limited diffusion radiative transfer.
In the stationary case this yields a Helmholtz equation

.. math::

   \kappa\Psi - \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi\right) = \frac{\eta}{c},

which is solved by a geometric multigrid method. The default boundary conditions are of the Robin type.
For fully transient radiative transport, we offer discretizations based on the backward Euler and TGA schemes. 

Monte Carlo methods
-------------------

All types of moment-closed radiative transfer equations contain nonphysical artifacts (which may or may not be acceptable).
For example, in the diffusion approximation the radiative flux is :math:`F = -\frac{c}{3\kappa}\nabla \Psi`, implying that photons can leak around boundaries.
I.e. the diffusion approximation does not correctly describe shadows.
It is possible to go beyond the diffusion approximation by also solving for higher-order moments like the radiative flux.
While such methods can describe shadows, they contain other nonphysical features.

.. figure:: figures/rte_comp.png
   :width: 720px
   :align: center

   Qualitative comparison between predictions made with a diffusion RTE solver and a Monte Carlo RTE solver.
   Left: Source term.
   Middle: Solution computed in the diffusion approximation with homogeneous Robin boundary conditions.
   Right: Solution computed with a Monte Carlo method. 

Both ''stationary'' and transient Monte Carlo methods are offered as an alternative to the diffusion approximation. 

photon particle
_______________

The Îto particle is a computational particle class in `PlasmaC` which can be used together with the particle tools in `Chombo`.
The following data fields are implemented in the particle:

.. code-block:: c++
   
   RealVect m_position;
   RealVect m_velocity;
   Real m_mass;
   Real m_kappa;

To obtain the fields, the user will call

.. code-block:: c++

   RealVect& position();
   RealVect& velocity();
   Real& mass();
   Real& diffusion();


All functions also have ``const`` versions.
Note that the field ``m_mass`` is the same as the *weight* of the computational particle.
The following functions are used to set the various properties:

.. code-block:: c++

   setPosition(const RealVect a_pos);
   setVelocity(const RealVect a_vel);
   setMass(const Real a_mass);
   setDiffusion(const Real a_diffusion;

Interaction with boundaries
___________________________


Stationary Monte Carlo
______________________

The stationary Monte Carlo method proceeds as follows.

1. For each cell in the mesh, draw a discrete number of photons :math:`\mathcal{P}\left(\eta \Delta V\Delta t\right)` where :math:`\mathcal{P}` is a Poisson distribution. The user may also choose to use pseudophotons rather than physical photons by modifying photon weights. Each photon is generated in the cell centroid :math:`\mathbf{x}_0` and given a random propagation direction :math:`\mathbf{n}`.

2. Draw a propagation distance :math:`r` by drawing random numbers from an exponential distribution :math:`p(r) = \kappa \exp\left(-\kappa r\right)`. The absorbed position of the photon is :math:`\mathbf{x} = \mathbf{x}_0 + r\mathbf{n}`.

3. Check if the path from :math:`\mathbf{x}_0` to :math:`\mathbf{x}` intersects an internal or domain boundary. If it does, absorb the photon on the boundary. If not, move the photon to :math:`\mathbf{x}` or reflect it off symmetry boundaries. 

4. Rebin the absorbed photons onto the AMR grid. This involves parallel communication. 

5. Compute the resulting photoionization profile. The user may choose between several different deposition schemes (like e.g. cloud-in-cell).
      
The Monte Carlo methods use computational particles for advancing the photons in exactly the same way a Particle-In-Cell method would use them for advancing electrons. Although a computational photon would normally live on the finest grid level that overlaps its position, this is not practical for all particle deposition kernels. For example, for cloud-in-cell deposition schemes it is useful to have the restrict the interpolation kernels to the grid level where the particle lives. In Chombo-speak, we therefore use a buffer region that extends some cells from a refinement boundary where the photons are not allowed to live. Instead, photons in that buffer region are transferred to a coarser level, and their deposition clouds are first interpolated to the fine level before deposition on the fine level happens. Selecting a deposition scheme and adjusting the buffer region is done through an input script associated with the solver. 
   
Transient Monte Carlo
_____________________

The transient Monte Carlo method is almost identical to the stationary method, except that it does not deposit all generated photons on the mesh but tracks them through time. The transient method is implemented as follows:

1. For each cell in the mesh, draw a discrete number of photons :math:`\mathcal{P}\left(\eta \Delta V\Delta t\right)` as above, and append these to the already existing photons. Each photon is given a uniformly distributed random creation time within :math:`\Delta t`. 
   
2. Each photon is advanced over the time step :math:`\Delta t` by a sequence of :math:`N` substeps (:math:`N` may be different for each photon).

   a. We compute :math:`N` such that we sample :math:`N\Delta \tau = \Delta t` with :math:`c\kappa\Delta\tau < 1`.

   b. A photon at position :math:`\mathbf{x}_0` is moved a distance :math:`\Delta \mathbf{x} = c\mathbf{n}\Delta\tau`. For each step we compute the absorption probability :math:`p = \kappa\left|\Delta\mathbf{x}\right|` where :math:`p\in[0,1]` is a uniform random number. If the photon is absorbed on this interval, draw a new uniform random number :math:`r \in [0,1]` and absorb the photon at the position :math:`\mathbf{x}_0 + r\Delta\mathbf{x}`. If the photon is not absorbed, it is moved to position :math:`\mathbf{x}_0 + r\Delta\mathbf{x}`.

3. Check if the path from :math:`\mathbf{x}_0` to :math:`\mathbf{x}` intersects an internal or domain boundary. If it does, absorb the photon on the boundary. If not, move the photon to :math:`\mathbf{x}`.

4. Rebin the absorbed photons onto the AMR grid. This involves parallel communication. 

5. Compute the resulting photoionization profile. The user may choose between several different deposition schemes (like e.g. cloud-in-cell).

Limitations
-----------
