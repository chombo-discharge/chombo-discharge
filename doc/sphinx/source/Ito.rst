.. _Chap:ItoDiffusion:

Îto diffusion
=============

The Ito diffusion model advances computational particles as Brownian walkers with drift:

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

The Îto particle is a computational particle class in `PlasmaC` which can be used together with the particle tools in `Chombo`.
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

Limitations
-----------
