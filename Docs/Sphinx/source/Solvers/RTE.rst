.. _Chap:RadiativeTransfer:

Radiative transfer
==================

.. _Chap:RtSolver:   

RtSolver
--------

Radiative transfer solvers are supported in the form of

* Diffusion solvers, i.e. first order Eddington solvers, which take the form of a Helmholtz equation.
* Particle solvers, which track photons are particles (e.g., Monte Carlo sampled solvers).
  
The solvers share a parent class ``RtSolver``, and code that uses only the ``RtSolver`` interface will should be able to switch between the two implementations.
Note, however, that the radiative transfer equation is inherently deterministic while Monte Carlo photon transport is inherently stochastic. 
The diffusion approximation relies on solving an elliptic equation in the stationary case and a parabolic equation in the time-dependent case, while the Monte-Carlo approach solves solves for fully transient or instantaneous transport.

.. tip::
   
   The source code for the solver is located in :file:`$DISCHARGE_HOME/Source/RadiativeTransfer` and it is a fairly lightweight abstract class.
   As with other solvers, ``RtSolver`` can use a specified :ref:`Chap:Realm`.

To use the ``RtSolver`` interface the user must cast from one of the inherited classes (see :ref:`Chap:DiffusionRTE` or :ref:`Chap:MonteCarloRTE`).
Since most of the ``RtSolver`` is an interface which is implemented by other radiative transfer solvers, documentation of boundary conditions, kernels and so on are found in the implementation classes.

.. _Chap:RtSpecies:

RtSpecies
_________

The class ``RtSpecies`` is an abstract base class for parsing necessary information into radiative transfer solvers.
When creating a radiative transfer solver one will need to pass in a reference to an ``RtSpecies`` instantiation such that the solvers can look up the required information.
Currently, ``RtSpecies`` is a lightweight class where the user needs to implement the function

.. literalinclude:: ../../../../Source/RadiativeTransfer/CD_RtSpecies.H
   :language: c++
   :lines: 48-53
   :dedent: 2

This absorption coefficient is used in both the diffusion (see :ref:`Chap:DiffusionRTE`) and Monte Carlo (see :ref:`Chap:MonteCarloRTE`) solvers.

.. important::
   
   Upon construction, one must set the class member ``m_name``, which is the name passed to the actual solver.

Setting the source term
_______________________

``RtSolver`` stores a source term :math:`\eta` on the mesh, which describes the number of photons that are generated produced per unit volume and time.
This variable can be set through the following functions:

.. literalinclude:: ../../../../Source/RadiativeTransfer/CD_RtSolver.H
   :language: c++
   :lines: 254-273
   :dedent: 2

The usage of :math:`\eta` varies between the different solvers.
It is possible, for example, to generate computational photons (particles) using :math:`\eta` when using Monte Carlo sampling, but this is not a requirement.

.. _Chap:DiffusionRTE:

Diffusion approximation
-----------------------

.. _Chap:EddingtonSP1:

EddingtonSP1
____________

The first-order diffusion approximation to the radiative transfer equation is encapsulated by the ``EddingtonSP1`` class which implements a first order Eddington approximation of the radiative transfer equation.
``EddingtonSP1`` implements ``RtSolver`` using both stationary and transient advance methods (e.g. for stationary or time-dependent radiative transport).
The source code is located in :file:`$DISCHARGE_HOME/RadiativeTransfer`. 

Equation(s) of motion
_____________________

In the diffusion approximation, the radiative transport equation is

.. math::
   :label: TransientDiffusionRTE

   \partial_t\Psi + \kappa\Psi - \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi\right) = \frac{\eta}{c},

where :math:`\Psi` is the radiative intensity (i.e., photons absorbed per unit volume`.
Here, :math:`\kappa` is the absorption coefficient (i.e., inverse absorption length).
This value can be spatially dependent, and is passed in through the :ref:`Chap:RtSpecies` function ``getAbsorptionCoefficient`` that was discussed above.
Note that in the context below, :math:`\kappa` is *not* the volume fraction of a grid cell but the absorption coefficient.
The above equation is called the Eddington approximation, with the closure relation being that the radiative flux is given by :math:`F = -\frac{c}{3\kappa}\nabla \Psi`.

In the stationary case this reduces to a Helmholtz equation

.. math::
   :label: StationaryDiffusionRTE

   \kappa\Psi - \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi\right) = \frac{\eta}{c},

Implementation
______________

``EddingtonSP1`` uses multigrid methods for solving :eq:`TransientDiffusionRTE` and :eq:`StationaryDiffusionRTE`, see :ref:`Chap:LinearSolvers`.
To advance the solution, one will call the member function

.. literalinclude:: ../../../../Source/RadiativeTransfer/CD_EddingtonSP1.H
   :language: c++
   :lines: 80-89
   :dedent: 2

Internally, this version will perform one of the following:

#. Solve :eq:`TransientDiffusionRTE` if using a *transient* solver.
   This is done using a backward Euler solve:

   .. math::

      \left(1+ \kappa \Delta t\right)\Psi^{k+1} - \Delta t \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi^{k+1}\right) = \Psi^{k} + \frac{\Delta t\eta^{k+1}}{c},

   This equation is a Helmholtz equation for :math:`\Psi^{k+1}` which is solved using geometric multigrid, see :ref:`Chap:LinearSolvers`.

#. Solve :eq:`StationaryDiffusionRTE` if using instantaneous photon transport.
   This is done directly with a geometric multigrid solver, see :ref:`Chap:LinearSolvers`.

.. _Chap:EddingtonSP1BC:
   
Boundary conditions
___________________

Simplified domain boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to set the following *simplified* boundary conditions on domain faces and embedded boundaries:
All of these boundary condition specifications take the form ``<type> <value>``.

#. Dirichlet, with a fixed value of :math:`\Phi`. E.g., ``dirichlet 0.0``.
#. Neumann, using a fixed value of :math:`\partial_n\Phi`. E.g., ``neumann 0.0``.
#. A *Larsen-type* boundary condition, which is an absorbing boundary condition in the form

   .. math::

      \kappa\partial_n\Psi + \frac{3\kappa^2}{2}\frac{1-3r_2}{1-2r_1}\Psi = g,

   where :math:`r_1` and :math:`r_2` are reflection coefficients and :math:`g` is a surface source, see :cite:`Larsen2002`.
   Note that when the user specifies the boundary condition value (e.g. by setting the BC function), he is setting the surface sourge :math:`g`.
   In the majority of cases, however, we will have :math:`r_1 = r_2 = g = 0` and the BC becomes

   .. math::

      \partial_n\Psi + \frac{3\kappa}{2}\Psi = g.   
      
   The user must then pass a value ``larsen <value>``, where the ``value`` corresponds to the souce term :math:`g`.
   Typically, this term is zero.

.. tip::
   
   For radiative transfer, the Larsen boundary condition is usually the correct one as it approximately describes outflow of photons on the boundary.
   In this case the correct boundary condition is ``larsen 0.0``.

Custom domain boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to use more complex boundary conditions by passing in ``dirichlet_custom``, ``neumann_custom``, or ``larsen_custom`` options through the solver configuration options (see :ref:`Chap:EddingtonInputOptions`).
In this case the ``EddingtonSP1`` solver will use a specified function at the domain edge/face which can vary spatially (and with time).
To specify that function, ``EddingtonSP1`` has a member function

.. literalinclude:: ../../../../Source/RadiativeTransfer/CD_EddingtonSP1.H
   :language: c++
   :lines: 115-125
   :dedent: 2

Here, the ``a_function`` argument is simply an alias:

.. literalinclude:: ../../../../Source/RadiativeTransfer/CD_EddingtonSP1DomainBc.H
   :language: c++
   :lines: 43-46
   :dedent: 2

Note that the boundary condition *type* is still Dirichlet, Neumann, or Larsen (depending on whether or not ``dirichlet_custom``, ``neumann_custom``, or ``larsen_custom`` was passed in). 
For example, to set the boundary condition on the left :math:`x` face in the domain, one can create a ``EddingtonSP1DomainBc::BcFunction`` object as follows:

.. code-block:: c++

   // Assume this has been instantiated. 
   RefCountedPtr<EddingtonSP1> eddingtonSolver;

   // Make a lambda which we can bind to std::function. 
   auto myValue = [](const RealVect a_pos, const Real a_time) -> Real {
      return a_pos[0] * a_time;
   }

   // Set the domain bc function in the solver. 
   eddingtonSolver.setDomainSideBcFunction(0, Side::Lo, myValue);

.. warning::

   A run-time error will occur if the user specifies one of the custom boundary conditions but does not actually set the function.
   
Embedded boundaries
^^^^^^^^^^^^^^^^^^^

On the EB, we currently only support constant-value boundary conditions.
In the input script, the user can specify

* ``dirichlet <value>`` For setting a constant Dirichlet boundary condition everywhere. 
* ``neumann <value>`` For setting a constant Neumann boundary condition everywhere. 
* ``larsen <value>`` For setting a constant Larsen boundary condition everywhere. 

The specification of these boundary conditions occurs in precise analogy with the domain boundary conditions, and are therefore not discussed further here.

.. _Chap:EddingtonInputOptions:

Solver configuration
____________________

The ``EddingtonSP1`` implementation has a number of configurable options for running the solver, and these are given below:

.. literalinclude:: ../../../../Source/RadiativeTransfer/CD_EddingtonSP1.options
   :emphasize-lines: 4, 6-9, 20-33
   :caption: ``EddingtonSP1`` solver configuration options. Run-time configurable options are highlighted.

The multigrid options are analogous to the multigrid options for :ref:`Chap:FieldSolverGMG`, see :ref:`Chap:MultigridTuning`.

.. _Chap:MonteCarloRTE:

Monte Carlo solver
------------------

``McPhoto`` defines a class which can solve radiative transfer problems using discrete photons.
The class derives from :ref:`Chap:RtSolver` and can thus be used also be used by applications that only require the :ref:`Chap:RtSolver` interface.
``McPhoto`` can provide a rather complex interaction with boundaries, such as computing the intersection between a photon path and a geometry, and thus capture shadows (which :ref:`Chap:EddingtonSP1` can not).

The Monte Carlo sampling is a particle-based radiative transfer solver, and particle-mesh operations (see :ref:`Chap:ParticleMesh`) are thus required in order to deposit the photons on a mesh if one wants to compute mesh-based absorption profiles.

.. tip::

   The ``McPhoto`` class is defined in :file:`$DISCHARGE_HOME/Source/RadiativeTransfer/CD_McPhoto.H`.
   See the `McPhoto C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classMcPhoto.html>`_ for further details.

The solver has multiple data holders for systemizing photons, which is especially useful during transport kernels where some of the photons might strike a boundary:

* In-flight photons.
* Bulk-absorbed photons, i.e., photons that were absorbed on the mesh.
* EB-absorbed photons, i.e., photons that struck the EB during a transport step.
* Domain-absorbed photons, i.e., photons that struck the domain edge/face during a transport step.
* Source photons, for letting the user pass in externally generated photons into the solver.

Various functions are in place for obtaining these particles:

.. literalinclude:: ../../../../Source/RadiativeTransfer/CD_McPhoto.H
   :language: c++
   :lines: 398-431
   :dedent: 2

Photon particle
_______________

The ``Photon`` particle is a simple encapsulation of a computational photon which is used by ``McPhoto``.
It derives from ``GenericParticle<2,1>`` and stores (in addition to the particle position):

* The particle weight.
* The particle mean absorption coefficient.
* The particle velocity/direction.

.. tip::

   The ``Photon`` class is defined in :file:`$DISCHARGE_HOME/Source/RadiativeTransfer/CD_Photon.H`

When defining the ``McPhoto`` class, the particle's absorption coefficient can be computed from the implementation of the absorption function method in :ref:`Chap:RtSpecies`.

Generating photons
__________________

There are several ways users can generate computational photons that are to be transported by the solver.

#. Fetch the *source photons* by calling

   .. literalinclude:: ../../../../Source/RadiativeTransfer/CD_McPhoto.H
      :language: c++
      :lines: 426-431
      :dedent: 2

   The source photons can then be filled and added to the other photons.

#. Add photons directly, by first obtaining the in-flight photons through

   .. literalinclude:: ../../../../Source/RadiativeTransfer/CD_McPhoto.H
      :language: c++
      :lines: 398-403
      :dedent: 2

   Photons can then be added directly.

#. If the source term :math:`\eta` has been filled, the user can call ``McPhoto::advance`` to have the solver generate the computational photons and then advance them.
   This is the correct approach for, e.g., applications that always use mesh-based photon source terms and want to have the computational photons be generated on the fly.
   
   .. warning::

      The ``advance`` function is *only* meant to be used together with a mesh-based source term that the user has filled prior to calling the method.

      When using the ``advance``, the number of photons that are generated are limit to a user-specified number (see :ref:`Chap:McPhotoOptions` for further details).

Transport modes
_______________

``McPhoto`` can be run as a fully transient, in which photons are tracked in time, or as an instantaneous solver.
For the instantaneous mode, photon absorption positions are stochastically sampled with Monte Carlo procedure and the photons are immediately absorbed on the mesh.
For the transient mode the photon advancement occurs over :math:`\Delta t`, so there is a limited distance (:math:`c \Delta t`) that the photons can propagate.
In this case, only some of the photons will be absorbed on the mesh whereas the rest may continue their propagation.

Instantaneous transport
^^^^^^^^^^^^^^^^^^^^^^^

When using instantaneous transport, any photon generated in a time step is immediately absorbed on the boundary through the following steps:

#. Optionally, have the solver generate photons to be transport (or add them externally).
#. Draw a propagation distance :math:`r` by drawing random numbers from an exponential distribution :math:`p(r) = \kappa \exp\left(-\kappa r\right)`.
   Here, :math:`\kappa` is computed by calling the underlying :ref:`Chap:RtSpecies` absorption function.
   The absorbed position of the photon is set to :math:`\mathbf{x} = \mathbf{x}_0 + r\mathbf{n}`.

   .. warning::

      In instantaneous mode photons might travel infinitely long, i.e. there is no guarantee that :math:`c\Delta t \leq r`.
#. Deposit the photons on the mesh.

Transient transport
^^^^^^^^^^^^^^^^^^^

The transient Monte Carlo method is almost identical to the stationary method, except that it does not deposit all generated photons on the mesh but tracks them through time.
For each photon, do the following:

#. Compute an absorption length :math:`r` by sampling the absorption function at the current photon position.
   
#. Each photon is advanced over the time step :math:`\Delta t` such that the position is

   .. math::

      \mathbf{x} = \mathbf{x}_0 + \mathbf{c}\Delta t.

#. Check if :math:`\left|\mathbf{x}-\mathbf{x}_0\right| < r` and if it is, absorb the photon on the mesh.
   
Other transport kernels
^^^^^^^^^^^^^^^^^^^^^^^

In addition to the above two methods, the solver interface permits users to add e.g. source photons externally and add them to the solvers' transport kernel. 

.. _Chap:McPhotoOptions:

Solver configuration
____________________


``McPhoto`` can be configured through its input options, see below:

.. literalinclude:: ../../../../Source/RadiativeTransfer/CD_McPhoto.options

.. tip::

   The ``McPhoto`` class includes a hidden input parameter ``McPhoto.dirty_sampling = true/false`` which enables a cheaper sampling method for discrete photons when calling the ``advance`` method.
   The caveat is that the method does not incorporate boundary intersect, only works for instantaneous propagation, and avoids filling the data holders that are necessary for load balancing.

Clarifications
^^^^^^^^^^^^^^

When computational photons are generated through the solver, users might have filled the source term differently depending on the application.
For example, users might have filled the source term with the number of photons generated per unit volume and time, or the *physical* number of photons to be generated. 
The two input options ``McPhoto.photon_generation`` and ``McPhoto.source_type`` contain the necessary specifications for ensuring that the user-filled source term can be translated properly for ensuring that the correct number of physical photons are generated.
Firstly, ``McPhoto.source_type`` contains the specification of what the source term contains:

* ``number`` if the source term contains the physical number of photons.
* ``volume`` if the source terms contains the physical number of photons generated per unit volume.
* ``volume_rate`` if the source terms contains the physical number of photons generated per unit volume and time.
* ``rate`` if the source terms contains the physical number of photons generated per unit time.

When ``McPhoto`` calculates the number of physical photons in a cell, it will automatically determine from ``McPhoto.source_type``, :math:`\Delta V` and :math:`\Delta t` how many physical photons are to be generated in each grid cell.

``McPhoto.photon_generation`` permits the user to turn on/off Poisson sampling when determining how many photons will be generated.
If this is set to *stochastic*, the solver will first compute the number of physical photons :math:`\overline{N}_\gamma^{\text{phys}}` following the procedure above, and then run a Poisson sampling such that the final number of physical photons is

.. math::

   N_{\gamma}^{\text{phys}} = P\left(\overline{N}_{\gamma}^{\text{phys}}\right).

Otherwise, if ``McPhoto.photon_generation`` is set to *deterministic* then the solver will generate

.. math::

   N_{\gamma}^{\text{phys}} = \overline{N}_{\gamma}^{\text{phys}}

photons.
Again, these elements are important because users might have chosen to perform the Poisson sampling outside of ``McPhoto``.

.. important::

   All of the above procedures are done *per-cell*.
  

Example application
-------------------

Example applications that use ``RtSolver`` are found in:

* :ref:`Chap:RadiativeTransferModel`.
* :ref:`Chap:CdrPlasmaModel`.
* :ref:`Chap:ItoKMC`.    
