.. _Chap:Electrostatics:
   
Electrostatic solver
====================

Here, we discuss the discretization of the equation

.. math::
   :label: Poisson

   \nabla\cdot\left(\epsilon_r\nabla\Phi\right) = -\frac{\rho}{\epsilon_0}


where :math:`\Phi` is the electric potential, :math:`\rho` is the space charge density, and :math:`\epsilon_0` is the vacuum permittivity.
The relative permittivity is :math:`\epsilon_r = \epsilon_r\left(\mathbf{x}\right)` and can additionally be discontinuous at gas-dielectric interfaces.

.. note::
   
   All current electrostatic field solvers solve for the potential at the cell center (not the cell centroid).
   The code for the electrostatics solver is given in :file:`/Source/Electrostatics` and :file:`/Source/Elliptic`.

.. _Chap:FieldSolver:

FieldSolver
-----------

``FieldSolver`` is an abstract class for electrostatic solves in an EB context and contains most routines required for setting up and solving electrostatic problems.
``FieldSolver`` can solve over three phases, gas, dielectric, and electrode, and thus it is uses ``MFAMRCellData`` functionality where data is defined over multiple phases (see :ref:`Chap:MeshData`).

Note that in order to separate the electrostatic solver interface from the implementation, ``FieldSolver`` is a pure class without knowledge of numerical discretizations.
Currently, our only supported subclass is :ref:`Chap:FieldSolverGMG`, which supplies a finite-volume discretization that is solved with geometric multigrid.

.. tip::
   
   See the `FieldSolver C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classFieldSolver.html>`_ for the complete interface.

On gas-dielectric interfaces we enforce an extra equation

.. math::
   :label: GaussBC
   
   \epsilon_1\partial_{n_1}\Phi + \epsilon_2\partial_{n_2}\Phi = \sigma/\epsilon_0

where :math:`\mathbf{n}_1 = -\mathbf{n}_2` are the normal vectors pointing away from interface, and :math:`\sigma` is the surface charge density.

We point out that this equation can be enforced in various formats. 
The most common case is that :math:`\partial_n\Phi` are free parameters and :math:`\sigma` is a fixed parameter.
However, we *can* also fix :math:`\partial_n\Phi` on one side of the boundary and let :math:`\sigma` be the free parameter.
When using ``FieldSolverGMG`` (see :ref:`Chap:FieldSolverGMG`), users can choose between these two natural boundary conditions, see :ref:`Chap:PoissonEBBC`.

Using FieldSolver
-----------------

Using the ``FieldSolver`` is usually straightforward by first constructing the solver and then parsing the class options.
This usually involves several steps, such as instantiating the solver and providing proper references to :ref:`Chap:AmrMesh` and the computational geometry.
In addition, the user must set up a voltage curve and associated boundary conditions. 

.. tip::
   
   A complete example of setting up electrostatic problems and solving them is given in :ref:`Chap:ElectrostaticsModel`.


Setting the voltage
-------------------

The function signature for setting the voltage on various parts is

.. literalinclude:: ../../../../Source/Electrostatics/CD_FieldSolver.H
   :language: c++
   :lines: 381-387
   :dedent: 2

This allows setting a time-dependent voltage on electrodes and domain boundaries.
As shown above, one can easily use ``std::function<Real(const Real)>`` or lambdas to set the voltage.
It is relatively straightforward to turn tabulated data (see :ref:`Chap:LookupTable`) into lambda functions.
A simple example that returns a linearly varying voltage is given below:

.. code-block:: c++

   FieldSolver* fieldSolver;
   
   Real myVoltage = [] (const Real a_time) -> Real {
      return 1.0*a_time;
   };

   fieldSolver->setVoltage(myVoltage);



.. _Chap:PoissonDomainBC:

Domain boundary conditions
--------------------------

Domain boundary conditions for the solver must be set by the user through an input script, whereas the boundary conditions on internal surfaces are Dirichlet by default.
Note that on multifluid-boundaries the boundary condition is enforced by the conventional matching boundary condition that follows from Gauss` law.

General format
______________

The most general form of setting domain boundary conditions for ``FieldSolver`` is to specify a boundary condition *type* (e.g., Dirichlet) together with a function specifying the value.
Domain boundary condition *types* are parsed through a member function ``FieldSolver::parseDomainBc``.
This function will read string identifiers from the input script, and these identifiers are either in the format ``<string> <float>`` (simplified format) or in the format ``<string>`` (general format). 
For setting general types of Neumann or Dirichlet BCs on the domain sides, one will specify 

.. code-block:: text

   FieldSolverGMG.bc.x.low  = dirichlet_custom
   FieldSolverGMG.bc.x.high = dirichlet_neumann   

Unfortunately, due to the many degrees of freedom in setting domain boundary conditions, the procedure is a bit convoluted.
We first explain the general procedure. 

``FieldSolver`` will always set individual space-time functions on each domain side, and these functions are always in the form

.. code-block:: c++

   std::function<Real(const RealVect a_position, const Real a_time)> bcFunction;

To set a domain boundary condition function on a side, one can use the following member function:

.. code-block:: c++

   void FieldSolver::setDomainSideBcFunction(const int a_dir,
		                             const Side::LoHiSide a_side,
					     const std::function<Real(const RealVect a_position, const Real a_time)>);

For a general way of setting the function value on the domain side, one will use the above function together with an identifier ``dirichlet_custom`` or ``neumann_custom`` in the input script.
This identifier simply tells ``FieldSolver`` to use that function to either specifiy :math:`\Phi` or :math:`\partial_n\Phi` on the boundary. 
These functions are then directly processed by the numerical discretizations.

.. note::

   On construction, ``FieldSolver`` will set all the domain boundary condition functions to a constant of one (because the functions need to be populated).

Simplified format
_________________

``FieldSolver`` also supports a simplified method of setting the domain boundary conditions, in which case the user will specify Neumann or Dirichlet values (rather than functions) for each domain side.
These values are usually, but not necessarily, constant values.

In this case one will use an identifier ``<string> <float>`` in the input script, like so:

.. code-block:: text

   FieldSolverGMG.bc.x.low  = neumann   0.0
   FieldSolverGMG.bc.x.high = dirichlet 1.0

The floating point number has a slightly different interpretation for the two types of BCs.
Moreover, when using the simplified format the function specified through ``setDomainSideBcFunction`` will be used as a multiplier rather than being parsed directly into the numerical discretization.
Although this may *seem* more involved, this procedure is usually easier to use when setting constant Neumann/Dirichlet values on the domain boundaries.
It also automatically provides a link between a specified voltage wave form and the boundary conditions (unlike the general format, where the user must supply that link themselves). 

Dirichlet
*********

When using simplified parsing of Dirichlet domain BCs, ``FieldSolver`` will generate and parse a different function into the discretizations.
This function is *not* the same function as that which is parsed through ``setDomainSideBcFunction``. 
In C++ pseudo-code, this function is in the format

.. code-block:: c++

   Real dirichletFraction;
   
   auto f = [&func, ...](const RealVect a_pos, const Real a_time) -> Real {
      return func(a_pos, a_time) * voltage(a_time) * dirichletFraction;
   };

where ``voltage`` is the voltage wave form specified through ``FieldSolver::setVoltage``, and ``dirichletFraction`` is a placeholder for the floating point number specified in the input script, i.e. the floating point number in the input option.
For Dirichlet boundary conditions the solver will always multiply the provided input function by the voltage waveform.
That is, the function ``func(a_pos, a_time)`` is the space-time function set through ``setDomainSideBcFunction``.
Recall that, by default, this function is set to one so that the default voltage that is parsed into the numerical discretization is simply the specified voltage multiplied by the specified fraction in the input script.
For example, using

.. code-block:: text

   FieldSolverGMG.bc.y.low  = dirichlet 0.0
   FieldSolverGMG.bc.y.high = dirichlet 1.0

will the set voltage on the lower y-plane to ground and the voltage on the upper y-plane to the live voltage.
Specifically, on the upper y-plane this specification will generate a potential boundary condition function of the type

.. code-block:: c++

   auto func = [](const RealVect a_pos, const Real a_time) {return 1.0};
   dirichletFraction = 1.0;
   
   auto bc = [func](const RealVect a_pos, const Real a_time) {
      return func(a_pos, a_time) * voltage(a_time) * dirichletFraction;
   };

In order to set the voltage on the domain side to also be spatially dependent, one can either use ``dirichlet_custom`` as an input option, or ``dirichlet <float>`` and set a different multiplier on the domain edge (face).
As an example. by specifying ``bc.y.high = dirichlet 1.234`` in the input script AND setting the multiplier on the wall as follows:

.. code-block:: c++

   auto wallFunc = [](const RealVect a_pos, const Real a_time) -> Real {
      return 1.0 - a_pos[1];
   };
   
   fieldSolver->setDomainSideBcFunction(1, Side::Hi, wallFunc);

we end up with a voltage of

.. math::

   V(\mathbf{x},t) = 1.234(1-y)V(t)

on the upper y-plane. 

Neumann
*******

When using simplified parsing of Neumann boundary conditions, the procedure is precisely like that for Dirichlet boundary conditions *except* that multiplication by the voltage wave form is not made.
I.e. the boundary condition function that is passed into the numerical discretization is

.. code-block:: c++
		
   Real neumannFraction;
   
   auto func = [&func, ...](const RealVect a_pos, const Real a_time) -> Real {
      return func(a_pos, a_time) * neumannFraction;
   };

Note that since ``func`` is initialized to one, the floating point number in the input option directly specifies the value of :math:`\partial_n\Phi`. 

.. _Chap:PoissonEBBC:
   
EB boundary conditions
----------------------

Electrodes
__________

For the current ``FieldSolver`` the natural BC at the EB is Dirichlet with a specified voltage, whereas on dielectrics we enforce :eq:`GaussBC`.
The voltage on the electrodes are automatically retrieved from the specified voltages on the electrodes in the geometry being used (see :ref:`Chap:ComputationalGeometry`).
The exception to this is that while :ref:`Chap:ComputationalGeometry` specifies that an electrode will be at some fraction of a specified voltage, ``FieldSolverGMG`` uses this fraction *and* the specified voltage wave form in ``setVoltage``.

To understand how the voltage on the electrode is being set, we first remark that our implementation uses a completely general specification of the voltage on each electrode in both space and time.
This voltage has the form

.. math::

   V_i = V_i\left(\mathbf{x}, t\right). 

where :math:`V_i` is the voltage on electrode :math:`i`.
It is possible to interact with this function directly, going through all electrodes and setting the electrode to be spatially and temporally varying.
The member function that does this is

.. literalinclude:: ../../../../Source/Electrostatics/CD_FieldSolver.H
   :language: c++
   :lines: 401-407
   :dedent: 2

Here, the type ``ElectrostaticEbBc::BcFunction`` is just an alias of ``std::function<Real(const RealVect a_position, const Real a_time)>``.
The voltage on an electrode :math:`i` could thus be set as

.. code-block:: c++

   int electrode;
   
   auto myElectrodeVoltage = [](const RealVect a_position, const Real a_time) -> Real{
       return 1.0;
   };

   fieldSolver->setElectrodeDirichletFunction(electrode, myElectrodeVoltage);

where the return value can be replaced by the user function.
In principle, one can then also set spatially varying voltages along an electrode.

In the majority of cases the voltage on electrodes is either a live voltage or ground.
Thus, although the above format is a general way of setting the voltage individually on each electrode (in both space and time) ``FieldSolver`` supports a simpler way of generating these voltage waveforms.
When ``FieldSolver`` is instantiated, it will interally generate these functions through simplified expression such that the user only needs to set a single wave form that applies to all electrodes.
The voltages that are set on the various electrodes are thus in the form:

.. code-block:: c++

   int electrode;
   Real voltageFraction;
   std::function<Real(const Real a_time)> voltageWaveForm;
   
   auto defaultElectrodeVoltage = [...](const RealVect a_position, const Real a_time) -> Real{
      return voltageFraction * voltageWaveForm(a_time); 
   };

   fieldSolver->setElectrodeDirichletFunction(electrode, defaultElectrodeVoltage);   

In summary, the default voltage which is set on an electrode is the voltage *fraction* specified on the electrodes (in :ref:`Chap:ComputationalGeometry`) multiplied by a voltage wave form (specified by ``FieldSolver::setVoltage``, as discussed above).


Dielectrics
___________

On dielectrics, we enforce the jump boundary condition directly, see :ref:`Chap:JumpCondition`.

Calling the solve function
--------------------------

The electrostatic solver in ``chombo-discharge`` has a lot of supporting functionality, but essentially relies on only one critical function:
Solving for the potential.
This is encapsulated by the pure member function

.. literalinclude:: ../../../../Source/Electrostatics/CD_FieldSolver.H
   :language: c++
   :lines: 107-117
   :dedent: 2

where ``a_phi`` is the resulting potential that was computing with the space charge density ``a_rho``, and surface charge density ``a_sigma``.
   
.. _Chap:FieldSolverGMG:   

FieldSolverGMG
--------------------

``FieldSolverGMG`` implements a multigrid routine for solving :eq:`Poisson`, and is currently the only implementation of ``FieldSolver``.

The discretization used by ``FieldSolverGMG`` is described in :ref:`Chap:LinearSolvers`.
The underlying solver type is a Helmholtz solver, but ``FieldSolverGMG`` considers only the Laplacian term.
For further details on the spatial discretization, see :ref:`Chap:LinearSolvers`.

Solver configuration
____________________

``FieldSolverGMG`` has a number of switches for determining how it operates.
Some of these switches are intended for parsing boundary conditions, whereas others are settings for operating multigrid or for I/O.
The current list of configuration options are indicated below

.. literalinclude:: ../../../../Source/Electrostatics/CD_FieldSolverGMG.options
   :language: text
   :emphasize-lines: 4,12-14,16-36
   :caption: Input options for the ``FieldSolverGMG`` class. Runtime adjustable options are highlighted.	      

Note that *all* options pertaining to IO or multigrid are run-time configurable (see :ref:`Chap:RuntimeConfig`).

Setting boundary conditions
___________________________

The flags that are in the format ``bc.coord.side`` (e.g., ``bc.x.low``) parse the domain boundary condition type to the solver.
See :ref:`Chap:PoissonDomainBC` for details.

The flag ``jump_bc`` indicates how the dielectric jump condition is enforced.
See :ref:`Chap:PoissonDielectricBC` for additional details.

.. note::
   
   Currently, we only solve the dielectric jump condition on gas-dielectric interfaces and dielectric-dielectric interfaces are not supported.
   If you want to use numerical mock-ups of dielectric-dielectric interfaces, you can change :math:`\epsilon_r` inside a dielectric, but note that the dielectric boundary condition :math:`\partial_{n_1}\Phi + \partial_{n_2}\Phi = \sigma/\epsilon_0` is *not* solved in this case.

Algorithmic adjustments
_______________________

By default, the Helmholtz operator uses a diagonally weighting of the operator using the volume fraction as weight.
This means that the quantity that is passed into ``AMRMultiGrid`` should be weighted by the volume fraction to avoid the small-cell problem of EB grids.
The flag ``kappa_source`` indicates whether or not we should multiply the right-hand side by the volume fraction before passing it into the solver routine.
If this flag is set to ``false``, it is an indication that the user has taken responsibility to perform this weighting prior to calling ``FieldSolver::solve(...)``.
If this flag is set to ``true``, ``FieldSolverGMG`` will perform the multiplication before the multigrid solve.

.. _Chap:MultigridTuning:

Tuning multigrid performance
____________________________

Multigrid operates by coarsening the solution (and the geometry with it) on a hierarchy of grid levels, and smoothing the solution on each level.
There are a number of factors that influence the multigrid performance.
Often the most critical factors are the radius of the cut-cell stencils and how far multigrid is allowed to coarsen.
In addition, the multigrid convergence is improved by increasing the number of smoothings per grid level (up to a certain point), as well as the type of smoother and bottom solver being used.
We explain these options below:

* ``FieldSolverGMG.gmg_verbosity``.
  Controls the multigrid verbosity.
  Setting it to a number :math:`> 0` will print multigrid convergence information.  
* ``FieldSolverGMG.gmg_use_default_settings``.
  Use default multigrid settings.
  This tends to make most problems converge.
* ``FieldSolverGMG.gmg_pre_smooth``.
  Controls the number of relaxations on each level during multigrid downsweeps.
* ``FieldSolverGMG.gmg_post_smooth``.
  Controls the number of relaxations on each level during multigrid upsweeps.
* ``FieldSolverGMG.gmg_bott_smooth``.
  Controls the number of relaxations before entering the bottom solve. 
* ``FieldSolverGMG.gmg_min_iter``.
  Sets the minimum number of iterations that multigrid will perform. 
* ``FieldSolverGMG.gmg_max_iter``.
  Sets the maximum number of iterations that multigrid will perform. 
* ``FieldSolverGMG.gmg_exit_tol``.
  Sets the exit tolerance for multigrid.
  Multigrid will exit the iterations if :math:`r < \lambda r_0` where :math:`\lambda` is the specified tolerance, :math:`r = |L\Phi -\rho|` is the residual and :math:`r_0` is the residual for :math:`\Phi = 0`.  
* ``FieldSolverGMG.gmg_exit_hang``.
  Sets the minimum permitted reduction in the convergence rate before exiting multigrid.
  Letting :math:`r^k` be the residual after :math:`k` multigrid cycles, multigrid will abort if the residual between levels is not reduce by at least a factor of :math:`r^{k+1} < (1-h)r^k`, where :math:`h` is the "hang" factor.
* ``FieldSolverGMG.gmg_min_cells``.
  Sets the minimum amount of cells along any coordinate direction for coarsened levels.
  Note that this will control how far multigrid will coarsen. Setting a number ``gmg_min_cells = 16`` will terminate multigrid coarsening when the domain has 16 cells in any of the coordinate direction. 
* ``FieldSolverGMG.gmg_bc_order``.
  Sets the stencil order for Dirichlet boundary conditions (on electrodes).
  Note that this is also the stencil radius. 
* ``FieldSolverGMG.gmg_bc_weight``. Sets the least squares stencil weighting factor for least squares gradient reconstruction on EBs.
  See :ref:`Chap:LeastSquares` for details. 
* ``FieldSolverGMG.gmg_jump_order``. Sets the stencil order when performing least squares gradient reconstruction on dielectric interfaces.
  Note that this is also the stencil radius. 
* ``FieldSolverGMG.gmg_jump_weight``.
  Sets the least squares stencil weighting factor for least squares gradient reconstruction on dielectric interfaces.
  See :ref:`Chap:LeastSquares` for details. 
* ``FieldSolverGMG.gmg_bottom_solver``.
  Sets the bottom solver type. 
* ``FieldSolverGMG.gmg_cycle``.
  Sets the multigrid method.
  Currently, only V-cycles are supported.
* ``FieldSolverGMG.gmg_smoother``.
  Sets the multigrid smoother.
* ``FieldSolverGMG.gmg_relax_factor``.
  Sets the multigrid relaxation factor.


  .. tip::
     
     Enabling this setting tends to make most problems converge quite well.


.. warning::

   When setting the bottom solver (which by default is a biconjugate gradient stabilized method) to a regular smoother, one must also specify the number of smoothings to perform.
   E.g., ``FieldSolverGMG.gmg_bottom_solver = simple 64``.
   Setting the bottom solver to ``simple`` without specifying the number of smoothings that will be performed will issue a run-time error. 
		

Adjusting output
________________

The user may plot the potential, the space charge, the electric, and the GMG residue as follows:

.. code-block:: text

   FieldSolverGMG.plt_vars  = phi rho E res     # Plot variables. Possible vars are 'phi', 'rho', 'E', 'res'

.. _Chap:PoissonDielectricBC:   

Saturation charge BC
____________________

As mentioned above, on dielectric interfaces the user can choose to specify which "form" of :eq:`GaussBC` to solve.
If the user wants the natural form in which the surface charge is the free parameter, he can specify

.. code-block:: text

   FieldSolverGMG.which_jump = natural

To use the other format (in which one of the fluxes is specified), use

.. code-block:: text

   FieldSolverGMG.which_jump = saturation_charge

.. note::
   
   The ``saturation_charge`` option will set the derivative of :math:`\partial_n\Phi` to zero on the gas side.
   Support for setting :math:`\partial_n\Phi` to a specified (e.g., non-zero) value on either side is missing, but is straightforward to implement.

Frequency dependent permittivity
--------------------------------

Frequency-dependent permittivities are fundamentally supported by the ``chombo-discharge`` elliptic discretization but none of the solvers implement it.
Recall that the polarization (in frequency space) is

.. math::

   \mathbf{P}(\omega) = \epsilon_0\chi(\omega)\mathbf{E}(\omega),

where :math:`\chi(\omega)` is the dielectric susceptibility.

There are two forms that ``chombo-discharge`` can support frequency dependent permittivities; through convolution or through auxiliary differential equations (ADEs).

Convolution approach
____________________

In the time domain, the displacement field is.

.. math::

   \mathbf{D}(t_k) = \epsilon_0\mathbf{E}(t_k) + \epsilon_0\int_0^{t_k} \chi(t)\mathbf{E}(t_k-t)\text{d}t.

There are various forms of discretizing the integral.
E.g. with the trapezoidal rule then

.. math::

   \begin{split}
   \int_0^{t_k}\chi(t)\mathbf{E}(t-t)\text{d}t &= \sum_{n=0}^{k-1} \int_{t_n}^{t_{n+1}}\chi(t)\mathbf{E}(t_k-t)\text{d}t \\
   &\approx \frac{1}{2}\sum_{n=0}^{k-1}\Delta t_n\left[\chi(t_n)\mathbf{E}(t_k-t_n) + \chi(t_{n+1})\mathbf{E}(t_k-t_{n+1})\right] \\
   &= \frac{\Delta t_0}{2}\chi_0\mathbf{E}(t_k) + \frac{1}{2}\sum_{n=1}^{k-1}\Delta t_n\chi_n\mathbf{E}(t_k-t_n) + \frac{1}{2}\sum_{n=0}^{k-1}\Delta t_n\chi_{n+1}\mathbf{E}(t_k-t_{n+1})
   \end{split}

The Gauss law becomes

.. math::

   \begin{split}
   \nabla\cdot\left[\left(1+\frac{\chi_0\Delta t_0}{2}\right)\mathbf{E}(t_k)\right] &= \frac{\rho(t_k)}{\epsilon_0}\\
   &- \nabla\cdot\left[\frac{1}{2}\sum_{n=1}^{k-1}\Delta t_n\chi_n\mathbf{E}(t_k-t_n) + \frac{1}{2}\sum_{n=0}^{k-1}\Delta t_n\chi_{n+1}\mathbf{E}(t_k-t_{n+1})\right].
   \end{split}

Note that the dispersion enters as an extra term on the right-hand side, emulating a space charge.
Unfortunately, inclusion of dispersion means that we must store :math:`\mathbf{E}(t_n)` for all previous time steps.

Auxiliary differential equation
_______________________________

With the ADE approach we seek a solution to :math:`\mathbf{P}(\omega) = \epsilon_0\chi(\omega)\mathbf{E}(\omega)` in the form

.. math::

   \sum_k a_k(i\omega)^k\mathbf{P}(\omega) = \epsilon_0\mathbf{E}(\omega),

where :math:`\sum a_k(i\omega)^k` is the Taylor series for :math:`1/\chi(\omega)`.
This can be written as a partial differential equation

.. math::
   
   \sum_{k}a_k\partial_t^k\mathbf{P}(t) = \epsilon_0\mathbf{E}(t).

This equation can be discretized using finite differences, and centering the solution on :math:`t_k` with backward differences yields an expression

.. math::
   
   \mathbf{P}^k = \epsilon_0C_0^k\mathbf{E}^k - \sum_{m>0} C_m^k\mathbf{P}^{k-m}.

where :math:`C_k` are stencil coefficients to be worked out for each case. 
The displacement field :math:`\mathbf{D}^k = \epsilon_0 \mathbf{E}^k + \mathbf{P}^k` is then

.. math::

   \mathbf{D} = \epsilon_0(1 + C_0^k)\mathbf{E} - \sum_{m>0} C_m^k\mathbf{P}^{k-m}.

The Gauss law yields

.. math::

   \nabla\cdot\left[\left(1 + C_0^k\right)\mathbf{E}^k\right] = \frac{\rho}{\epsilon_0} - \frac{1}{\epsilon_0}\nabla\cdot\sum_{m>0} C_m^k\mathbf{P}^{k-m}.

Unlike the convolution approach, this only requires storing terms required for the ADE description. 
This depends both on the order of the ADE, as well as it's discretization.
Normally, the ADE is a low-order PDE and a few terms are sufficient.    
   
Limitations
-----------

.. warning::

   There is currently a bug where having a dielectric interface align *completely* with a grid face will cause the cell to be identified as an electrode EB.
   This bug is due to the way ``Chombo`` handles cut-cells that align completely with a grid face.
   In this case the cell with volume fraction :math:`\kappa = 1` will be identified as an irregular cell.
   For the opposite phase (i.e., viewing the grids from inside the boundary) the situation is opposite and thus the two "matching cells" can appear in different grid patches.
   A fix for this is underway.
   In the meantime, a sufficient workaround is simply to displace the dielectric slightly away from the interface (any non-zero displacement will do).

Example application(s)
----------------------

Example applications that use the electrostatics capabilities are:

* :ref:`Chap:ElectrostaticsModel`.
* :ref:`Chap:CdrPlasmaModel`.   
