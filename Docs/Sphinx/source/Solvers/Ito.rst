.. _Chap:ItoDiffusion:

Îto diffusion
=============

The Îto diffusion model advances computational particles as Brownian walkers with drift:

.. math::
   d\mathbf{X}_i = \mathbf{v}_idt + \sqrt{2D_i}\mathbf{W}_i dt,

where :math:`\mathbf{X}_i` is the spatial position of a particle :math:`i`, :math:`\mathbf{v}_i` is the drift coefficient and :math:`D_i` is the diffusion coefficient *in the continuum limit*.
That is, both :math:`\mathbf{v}_i` and :math:`D_i` are the quantities that appear in :ref:`Chap:CDR`.
The vector term :math:`\mathbf{W}_i` is a Gaussian random field with a mean value of 0 and standard deviation of 1.

The code for Îto diffusion is given in :file:`/src/ito_solver` and only a brief explanation is given here.
The source code is used by a physics module in :file:`/physics/brownian_walker` and in the regression test :file:`/regression/brownian_walker`. 

.. _Chap:ito_particle:

The Îto particle
----------------

The Îto particle is a computational particle class in `chombo-discharge` which can be used together with the particle tools in `Chombo`.
The following data fields are implemented in the particle:

.. code-block:: c++
   
   RealVect m_position;
   RealVect m_velocity;
   Real m_mass;
   Real m_diffusion;

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

.. _Chap:ito_species:

ito_species
-----------

``ito_species`` is a class for parsing information into the solver class.
The constructor for the ``ito_species`` class is

.. code-block:: c++

   ito_species(const std::string a_name, const int a_charge, const bool a_mobile, const bool a_diffusive);

and this will set the name of the class, the charge, and whether or not the transport kernels account for drift and diffusion.

Setting initial conditions
__________________________

In order to set the initial conditions the user must fill the list ``List<ito_particle> m_initial_particles`` in ``ito_species``. 
When ``initial_data()`` is called from ``ito_solver``, the initial particles are transferred from the instance of ``ito_species`` and into the instance of ``ito_solver``.

We remark that it is a bad idea to replicate the initial particle list over all MPI ranks in a simulation.
If one has a list of initial particles, or wants to draw a specified number of particles from a distribution, the initial particles *must* be distributed over the available MPI ranks.
For example, the code in :file:`/physics/brownian_walker/brownian_walker_species.cpp` draws a specified number of particles distributed over all MPI ranks as (the code is called in ``brownian_walker_species::draw_initial_particles``)

.. code-block:: c++

  // To avoid that MPI ranks draw the same particle positions, increment the seed for each rank
  m_seed += procID();

  // Set up the RNG
  m_rng = std::mt19937_64(m_seed);
  m_gauss = std::normal_distribution<Real>(0.0, m_blob_radius);
  m_udist11 = std::uniform_real_distribution<Real>(-1., 1.);

  // Each MPI process draws the desired number of particles from a distribution
  const int quotient  = m_num_particles/numProc();
  const int remainder = m_num_particles % numProc();
  
  Vector<int> particlesPerRank(numProc(), quotient);
  
  for (int i = 0; i < remainder; i++){ 
    particlesPerRank[i] += 1;
  }

  // Now make the particles
  m_initial_particles.clear();
  for (int i = 0; i < particlesPerRank[procID()]; i++){
    const Real weight  = 1.0;
    const RealVect pos = m_blob_center + random_gaussian();
    m_initial_particles.add(ito_particle(weight, pos));
  }

Computing time steps
--------------------

The signatures for computing a time step for the ``ito_solver`` are given separately for the drift part and the diffusion part.

Drift
_____

The drift time step routines are implemented such that one restricts the time step such that the fastest particle does not move more than a specified number of grid cells. 

For the drift, the signatures are

.. code-block:: c++
		
  Real compute_min_drift_dt(const Real a_maxCellsToMove) const;
  Vector<Real> compute_drift_dt(const Real a_maxCellsToMove) const;
  
  Vector<Real> compute_drift_dt() const; // Compute dt on all AMR levels, return vector of time step
  Real compute_drift_dt(const int a_lvl) const;
  Real compute_drift_dt(const int a_lvl, const DataIndex& a_dit, const RealVect a_dx) const;

These last three functions all compute :math:`\Delta t = \Delta x/Max(v_x, v_y, v_z)` on the the various AMR levels and patches.
The routine

.. code-block:: c++

   Vector<Real> compute_drift_dt(const Real a_maxCellsToMove) const;

simply scales :math:`\Delta t` by ``a_maxCellsToMove`` on every level.
Finally, the function ``compute_min_drift_dt(...)`` computes the smallest time step across every AMR level. 

Diffusion
_________

The signatures for the diffusion time step are similar to the ones for drift:

.. code-block:: c++

   Real compute_min_diffusion_dt(const Real a_maxCellsToMove) const;
   Vector<Real> compute_diffusion_dt(const Real a_maxCellsToMove) const;

   Vector<Real> compute_diffusion_dt() const;
   Real compute_diffusion_dt(const int a_lvl) const;
   Real compute_diffusion_dt(const int a_lvl, const DataIndex& a_dit, const RealVect a_dx) const;

In these routines, the time step is computed as :math:`\Delta t = \frac{\Delta x}{\sqrt{2D}}`.
Note that there is still a chance that a particle jumps further than specified by ``a_maxCellsToMove`` since the diffusion hop is

.. math::

   \mathbf{d} = \sqrt{2D}\mathbf{Z}\Delta t,

where :math:`\mathbf{Z}` is a random Gaussian. 
The probability that a diffusion hop leads to a jump larger than :math:`N` cells can be evaluated and is :math:`P = \textrm{erf}\left(\sqrt{2}N\right)`. It is useful to keep this probability in mind when deciding on the PVR. 

Remapping particles
-------------------

Particle remapping has been implemented for the whole AMR hierarchy as a two step process.

1. Perform two-level remapping where particles are transferred up or down one grid level if they move out the level PVR.
2. Gather all particles that are remnant in the outcast list on the coarsest level, and then distribute them back to their appropriate levels. For example, particles that hopped over more than one refinement boundary cannot be transferred with a (clean) two-level transfer. 



Limitations
-----------

Example application
-------------------

An example application of usage of the ``ItoSolver`` is found in :ref:`Chap:BrownianWalkerModel`. 
