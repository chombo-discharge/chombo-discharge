.. _Chap:TemporalDiscretization:

Temporal discretization
=======================

In this chapter we discuss the supported temporal integrators for `PlasmaC`, and discuss their input parameters. These integrators differ in their level of efficiency and accuracy. Currently, none of the integrators can subcycle in time.

Time step limitations
---------------------

Our time integrators have different time step limitations. Fully explicit codes are limited by the advective and diffusive CFL constraints and usually also the dielectric relaxation time. However, some of the `PlasmaC` integrators eliminate the dielectric relaxation time, and all integrators can handle diffusion either implicitly or explicitly. In some cases only the advective CFL constraint is the only time step restriction.


  
Deterministic integrators
-------------------------

.. _Chap:godunov:

godunov
_______
The ``godunov`` temporal integrator is a rather unsophisticated, but very stable, temporal integrator. ``godunov`` uses an operator splitting between charge transport and plasma chemistry in the following way:

1. Advance :math:`\partial_t\phi = -\nabla\cdot\left(\mathbf{v}\phi - D\nabla\phi + \sqrt{2D\phi}\mathbf{Z}\right)` (and the surface charge solver) over a time step :math:`\Delta t`.
2. Compute the electric field 
3. Advance the plasma chemistry over the same time step using the field computed in 2). I.e. advance :math:`\partial_t\phi = S` over a time step :math:`\Delta t`.   
4. Move photons and deposit them on the mesh. 
      
Various integration options for the transport and chemistry steps are available but are discussed elsewhere. Note that the ``godunov`` integrator uses a semi-implicit coupling between the plasma chemistry terms and the electric field, and therefore eliminates the so-called dielectric relaxation time. The formal order of convergence of the ``godunov`` integrator is 1, but the accuracy can be quite good depending on the transport and chemistry schemes that are chosen.



.. _Chap:strang2:
    
strang2
_______
The ``strang2`` temporal integrator uses a second order Strang splitting between advection-reaction and diffusion, and supports up to fifth order Runge-Kutta schemes for the advection-reaction part. The ``strang2`` integrator integrates the equations of motion in the following way:

.. math::
   \phi^{k+1} = \exp\left(\frac{\Delta t}{2}\mathcal{F}_{\textrm{AR}}\right)\exp\left(\Delta t\mathcal{F}_{\textrm{D}}\right)\exp\left(\frac{\Delta t}{2}\mathcal{F}_{\textrm{AR}}\right)\phi^k,

where :math:`\mathcal{A}` is the advection-reaction operator and :math:`\mathcal{D}` is the diffusion operator. To perform the advection-reaction advance, we use of various forms strongly stability preserving Runge-Kutta (SSPRK) methods of order :math:`p` and stages :math:`s`. The SSPRK methods advance in the form

.. math::
   :nowrap:
   
   \begin{align}
   \phi^{(0)} &= \phi^k, \\
   \phi^{(i)} &= \sum_{k=0}^{i-1}\left[\alpha_{ik}\phi^{(k)} + \Delta t\mathcal{A}\left(\phi^{(k)}\right)\right], \quad i=1,2,\ldots, s,\\
   \phi^{(k+1)} &= \phi^{(s)}.
   \end{align}
      
In `PlasmaC`, charge injection is a part of the advective discretization so that the temporal discretization of :math:`\partial_t\sigma` follows the same discretization. 

The use of higher order methods for the advection-reaction advance allows use of longer time steps than for lower order methods for the same accuracy. However SSPRK methods with :math:`s=p` have a CFL constraint of 1, so that their relative efficiency is :math:`1/p` when compared to the explicit Euler method. Higher-order methods are therefore not always equivalent to better solver performance since the extra accuracy is traded for extra right hand side evaluations. We therefore explore the use of :math:`s>p` SSPRK methods that use more stages in order to gain efficiencies closer to the Euler method. The complete list of the SSPRK methods that we support, together with the CFL limits and effective CFL numbers are provided in the table below. For SSPRK methods of order 2 the CFL stability limit is given :math:`s-1` which yields an effective CFL of :math:`(s-1)/s`. For :math:`s=5`, for example, this gives a scheme that performs 60 percent faster than the equivalent two-stage method (provided one maximizes the CFL in each case). Likewise, the SSPRK(4,3) scheme has a CFL stability limit of 2, which gives an effective CFL of 0.5, being equivalent to the two-stage second-order method. However, the SSPRK(4,3) provides better accuracy of the advective advancement. 


==========  =============   ===============
Method      CFL             Effective CFL
==========  =============   ===============
SSPRK(s,2)  :math:`s-1`     :math:`(s-1)/s`
SSPRK(3,3)  :math:`1`       :math:`1/3` 
SSPRK(4,3)  :math:`2`       :math:`1/2` 
SSPRK(5,3)  :math:`2.65`    :math:`0.53`
SSPRK(5,4)  :math:`1.508`   :math:`0.3`
==========  =============   ===============

The diffusion advance uses the implicit scheme by Twizel, Gumel, and Arigu (TGA) and discretizes as

.. math::
   \left(I - \mu_1\mathcal{F}_{\textrm{D}}\right)\left(I - \mu_2\mathcal{F}_{\textrm{D}}\right)\phi^{k+1} = \left(I + \mu_3\mathcal{F}_{\textrm{D}}\right)\phi^k,

where :math:`\mu_1`, :math:`\mu_2`, :math:`\mu_3` are

.. math::
   :nowrap:
   
   \begin{align}
   \mu_1 &= \frac{a - \sqrt{a^2 - 4a + 2}}{2}\Delta t,\\
   \mu_2 &= \frac{a + \sqrt{a^2 - 4a + 2}}{2}\Delta t,\\
   \mu_3 &= (1-a)\Delta t. 
   \end{align}

Finally, note that the advective advance is performed with :math:`\Delta t/2`, and that the advancement can be made with twice the CFL number indicated in the table above. In order to estimate the numerical cost of the splitting method with :math:`s` Runge-Kutta stages, we remark that each Runge-Kutta stages requires one electric field update and one radiative transfer update for each RTE . In addition, there should be one electric field update after the diffusion update, and there will be two elliptic solves for each diffusive species. E.g. if only electrons are diffusive and we use a three-term RTE model , the :math:`s=4` method of order 3 will perform 35 elliptic updates per time step at a maximum CFL up to 4.

.. _Chap:SISDC:

IMEX SDC
________
``imex_sdc`` is a semi-implicit spectral deferred correction method for the `PlasmaC` equation set and is an adaptive high-order discretization with implicit diffusion. This method integrates the advection-diffusion-reaction equations in the following way.

Spectral deferred corrections
*****************************

Given an ordinary differential equation (ODE) as

.. math::
   \frac{\partial u}{\partial t} = F(u,t), \quad u(t_0) = u_0,

the exact solution is
  
.. math::
   u(t) = u_0 + \int_{t_0}^tF\left(u,\tau\right)d\tau.
   
Denote an approximation to this solution by :math:`\widetilde{u}(t)` and the correction by :math:`\delta(t) = u(t) - \widetilde{u}(t)`. The measure of error in :math:`\widetilde{u}(t)` is then

.. math::
   R(\widetilde{u}, t) = u_0 + \int_{t_0}^tF(\widetilde{u}, \tau)d\tau - \widetilde{u}(t).

Equivalently, since :math:`u = \widetilde{u} + \delta`, we can write

.. math::
   \widetilde{u} + \delta = u_0 + \int_{t_0}^t F\left(\widetilde{u}+\delta, \tau\right)d\tau. 

This yields

.. math::
   \delta = \int_{t_0}^t\left[F\left(\widetilde{u}+\delta, \tau\right) - F\left(\widetilde{u}, \tau\right)\right]d\tau + R\left(\widetilde{u},t\right). 

This is called the correction equation. The goal of SDC is to iteratively solve this equation in order to provide a high-order discretization. 

We now discuss the semi-implicit SDC (SISDC) method. First, we apply the method of lines (MOL) such that

.. math::
   :nowrap:
   
   \begin{align}
   \frac{d\phi_{\mathbf{i}}}{dt} &= \mathcal{F}_{\textrm{AR}}\left(t, \phi_{\mathbf{i}}\right) + \mathcal{F}_{\textrm{D}}\left(t, \phi_{\mathbf{i}}; \mathbf{E}_{\mathbf{i}}\right), \\
   \frac{d\sigma_{\mathbf{i}}}{dt} &= \mathcal{F}_{\sigma}\left(t, \phi_{\mathbf{i}}\right),
   \end{align}

where :math:`\phi_{\mathbf{i}}` denotes a cell-averaged variable, :math:`\mathcal{F}_{\sigma}` is as described in :ref:`Chap:SpatialDiscretization`, :math:`\mathcal{F}_{\textrm{AR}}\left(t, \phi_{\mathbf{i}}\right) = -D^c_{\mathbf{i}} + S_{\mathbf{i}}` is the advection-reaction operator , and :math:`\mathcal{F}_{\textrm{D}}(t, \phi_{\mathbf{i}}; \mathbf{E}_{\mathbf{i}}) = \frac{1}{\kappa_{\mathbf{i}}}\int_{V_{\mathbf{i}}}\left[\nabla\cdot\left(D\nabla\phi\right)\right]dV_{\mathbf{i}}` is the diffusion operator. Note that the advective operator contains the hybrid divergence discussed in :ref:`Chap:AdvectiveDiscretization` and :math:`\mathcal{F}_{\textrm{D}}` is parametrically coupled to :math:`\mathbf{E}` through :math:`D = D\left(\mathbf{E}\right)` (we use a semi-colon to indicate this dependence). Strictly speaking, :math:`\mathcal{F}_{\textrm{AR}}` is parametrically coupled in the same way through the  mobilities and boundary conditions, and additionally coupled to :math:`\Psi` through source terms so that the notation :math:`\mathcal{F}_{\textrm{AR}}\left(t, \phi_{\mathbf{i}}; \mathbf{E}_{\mathbf{i}}, \Psi_{\mathbf{i}}\right)` would be appropriate. However, charge injection, advection, and chemistry will be integrated explicitly so this dependence is notationally suppressed. On the other hand, the diffusion part will be solved with the backward Euler method - which yields a Helmholtz equation - and so we need to maintain this dependence for now. Later, we will clarify how this dependence is resolved. The rationale for solving diffusion implicitly is due to the numerical time step constraint of explicit diffusion methods which scales as :math:`\mathcal{O}\left(\Delta x^2\right)`, whereas advection scales more favorably at :math:`\mathcal{O}\left(\Delta x\right)`. We have chosen to integrate the reactive terms explicitly. The reason is that the reactive terms can be non-local, i.e. they can depend on the electron gradient. This is for example the case for fluid models in the local energy approximation where the electron energy source term contains terms that are proportional to the electron diffusion term :math:`D_e\nabla\phi_e`. Implicit discretization of the reactive terms then yield a fully coupled system rather than systems coupled only within individual cells. Charge injection is also handled explicitly. This design choice is mandated by the fact that implicit charge injection through the diffusion terms couples every diffusive species, leading to a system of diffusion equations that are fully coupled through their boundary conditions. Although charge injection could reasonably be performed as a separate step, this leads to numerical instabilities for cut-cell methods since the injected charge must be normalized by the volume fraction of the cell (which can be arbitrarily small). 

SISDC predictor
***************
Next, we present the SISDC method. In what follows, we suppress the index :math:`{\mathbf{i}}` as it is not explicitly needed. Given an interval :math:`[t_n, t_{n+1}]` on which a solution is sought, SDC methods divide this interval into :math:`p` subintervals :math:`t_n = t_{n,0} < t_{n,1} < \ldots < t_{n,p} = t_{n+1}`. Our discussion, however, pertains only to the interval :math:`[t_n, t_{n+1}]` so we compress the notation to :math:`t_m\equiv t_{n,m}`. We obtain an initial solution :math:`\phi_{m}^0, m=0,1,\ldots,p` as the semi-implicit advance

.. math::
   :nowrap:

   \begin{align}
   \phi_{m+1}^0 &= \phi_m^0 + \Delta t_m\left[\mathcal{F}_{\textrm{AR}}\left(t_m,\phi_m^0\right) + \mathcal{F}_{\textrm{D}}\left(t_{m+1},\phi_{m+1}^0; \mathbf{E}_{m+1}^0\right)\right],\\
   \sigma_{m+1}^0 &= \sigma_m^0 + \Delta t_mF_\sigma\left(t_m,\phi_m^0\right).
   \end{align}

This defines a Helmholtz problem for :math:`\phi_{m+1}^0` through :math:`\mathcal{F}_{\textrm{D}}`. Generally, the upper subscript denotes an SDC iteration where subscript 0 is the SISDC predictor, and we also have :math:`\phi_0^0 = \phi(t_n)` and :math:`\sigma_0^0 = \sigma(t_n)`. This predictor treats advection and chemistry terms explicitly, and diffusion implicitly. Other types of semi-implicit or multi-implicit couplings are possible :cite:`Bourlioux2003,Layton2004,Nonaka2012`. SDC improves this solution by using deferred corrections: Given a numerical solution :math:`\phi_{m+1}^k`, we compute an error :math:`\delta_{m+1}^k` and obtain the next iterate :math:`\phi_{m+1}^{k+1} = \phi_{m+1}^k + \delta_{m+1}^k`. Each iteration raises the discretization order by one :cite:`Dutt2000,Minion2003`, to maximum order :math:`p+1`. Critical to the success of this approach is the precise evaluation of the numerical quadrature. 

The parametric coupling of the electric field complicates things since the predictor contains :math:`\mathbf{E}_{m+1}^0 = \mathbf{E}\left(\phi_{m+1}^0\right)`, implying that the Poisson equation and the diffusion advance require concurrent solves for the diffusion update. We simplify this system by using a weak coupling by first computing

.. math::
   :nowrap:
   
   \begin{align}
   \phi_{m+1}^{0,\ast} &= \phi_m^0 + \Delta t_m\mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^0\right), \\
   \sigma_{m+1}^0 &= \sigma_m^0 + \Delta t_mF_\sigma\left(t_m, \phi_m^0\right),
   \end{align}

Next, we will approximate :math:`\mathbf{E}_{m+1}^{0}` for use in the predictor. There are two choices for this coupling; one may either use :math:`\mathbf{E}_m^0` for computation of the diffusion coefficients, which we will refer to as the semi-implicit coupling, or one may use fixed-point iteration and compute :math:`\mathbf{E}_{m+1}^{0,\ast} = \mathbf{E}\left(\phi_{m+1}^{0, \ast}, \sigma_{m+1}^0\right)`, followed by the diffusion advance

.. math::
   \phi_{m+1}^{0,\dagger} = \phi_{m+1}^{0,\ast} + \Delta t_m\mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^0; \mathbf{E}_{m+1}^\ast\right),

which we will refer to as the implicit coupling. This is e.g. the electric field coupling used in :cite:`Marskar2019`. This approximation can be improved by using more fixed-point iterations that computes :math:`\mathbf{E}_{m+1}^{0,\dagger} = \mathbf{E}\left(\phi_{m+1}^{0,\dagger}, \sigma_{m+1}^0\right)` and then re-solves the predictor equation with :math:`\mathbf{E}_{m+1}^{0,\dagger}` in place of :math:`\mathbf{E}_{m+1}^{0,\ast}`. The process can then be repeated for increased accuracy. Regardless of which coupling is used, we have now calculated :math:`\phi_{m+1}^0`, :math:`\sigma_{m+1}^0`, through which we obtain :math:`\mathbf{E}_{m+1}^0 = \mathbf{E}\left(\phi_{m+1}^0, \sigma_{m+1}^0\right)`, and :math:`\Psi_{m+1}^0 = \Psi\left(\mathbf{E}_{m+1}^0, \phi_{m+1}^0\right)`. Finally, we remark that the SISDC predictor is a sequentially advanced semi-implicit Euler method, which is locally second order accurate and globally first order accurate. Each step of the predictor can be thought of as a Godunov splitting between the advective-reactive and diffusive terms. 

SISDC corrector
***************
Next, the semi-implicit discretization of the correction equation is

.. math::
   \begin{split}
   \delta_{m+1}^k &= \delta_m^k  + \Delta t_m\left[\mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^k + \delta_m^k\right) - \mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^k\right)\right.\\
   &+ \left.\mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^k + \delta_{m+1}^k; \mathbf{E}_{m+1}^k\right) - \mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^k; \mathbf{E}_{m+1}^k\right)\right] - \left(R_{m+1}^k - R_{m}^k\right).
   \end{split}

We furthermore define

.. math::
   \begin{split}
   R_{m+1}^k - R_m^k &= \int_{t_m}^{t_{m+1}}\left[\mathcal{F}_{\textrm{AR}}\left(\phi^k\right) + \mathcal{F}_{\textrm{D}}\left(\phi^k; \mathbf{E}^k\right)\right]d\tau - \phi_{m+1}^k + \phi_m^k \\
   &\equiv I_m^{m+1}\left(\phi^k\right) - \phi_{m+1}^k + \phi_m^k. 
   \end{split} 

Evaluation of :math:`I_m^{m+1}` yields :math:`p` quadrature rules and we may write

.. math::
   I_m^{m+1}\left(\phi^k\right) = \sum_{l=0}^p q_m^l\left[\mathcal{F}_{\textrm{AR}}\left(t_l, \phi^k_l\right) + \mathcal{F}_{\textrm{D}}\left(t_l, \phi^k_l; \mathbf{E}_l^k\right)\right],

where the weights :math:`q_m^l` are quadrature weights. The final update for :math:`\phi^{k+1}_{m+1}` is then

.. math::
   \begin{split}
   \phi_{m+1}^{k+1} &= \phi_{m}^{k+1} + \Delta t_m\left[\mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^{k+1}\right) -\mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^{k}\right)\right.\\
   & + \left.\mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^{k+1}; \phi_{m+1}^{k+1}\right) - \mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^{k}; \mathbf{E}_{m+1}^k\right)\right] + I_{m}^{m+1}\left(\phi^k\right).
   \end{split}

With the exception of :math:`\mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^{k+1}; \mathbf{E}_{m+1}^{k+1}\right)`, all quantities on the right-hand are known and the correction equation is reduced to a Helmholtz equation for :math:`\phi_{m+1}^{k+1}` with error :math:`\delta_{m+1}^k = \phi_{m+1}^{k+1} - \phi_{m+1}^k`. An analogous equation is found for :math:`\sigma_{m+1}^{k+1}`.

The correction step has the same coupling to the electric field as the prediction step in that :math:`\mathbf{E}_{m+1}^{k+1}` appears in the update equation for :math:`\phi_{m+1}^{k+1}`. As for the prediction, we use a weak coupling through which we first compute

.. math::
   :nowrap:
   
   \begin{align}
   \phi_{m+1}^{k+1,\ast} &= \phi_m^{k+1} + \Delta t_m\left[\mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^{k+1}\right) - \mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^{k}\right)\right] + I_m^{m+1}\left(\phi^k\right),\\
   \sigma_{m+1}^{k+1} &= \sigma_m^{k+1} + \Delta t_m\left[F_\sigma\left(t_m, \phi_m^{k+1}\right) - F_\sigma\left(t_m, \phi_m^{k}\right)\right] + \Sigma_m^{m+1}\left(\phi^k\right). 
   \end{align}

The solution for :math:`\sigma_{m+1}^{k+1}` is final since all charge is injected through the advection operator for :math:`\phi`. The term :math:`\Sigma_m^{m+1}` contains the injected charge through :math:`I_m^{m+1}\left(\phi^k\right)`, as was discussed in :ref:`Chap:SpatialDiscretization`. We then solve

.. math::
   \phi_{m+1}^{k+1} = \phi_{m+1}^{k+1, \ast} + \Delta t_m\left[\mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^{k+1}; \mathbf{E}_{m+1}^{k+1}\right) - \mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^{k}; \mathbf{E}_{m+1}^k\right)\right],

with some approximation for :math:`\mathbf{E}_{m+1}^{k+1}`. As before, this coupling can be made either semi-implicitly or implicitly. The corrector equation defines a Helmholtz equation for :math:`\phi_{m+1}^{k+1}` using :math:`\phi_{m+1}^{k+1,\ast}` as the previous solution and :math:`-\mathcal{F}_{\textrm{D}}\left(\phi_{m+1}^{k}; \mathbf{E}_{m+1}^k\right)` as a source term.

Order, stability, and computational cost
****************************************
For consistency with the literature, denote the SISDC method which uses :math:`P` nodes (i.e. :math:`P-1` subintervals) and :math:`K` total iterations (i.e. :math:`K-1` iterations of the correction equation) by :math:`\verb|SISDC|_P^K`. This method will have a global order of accuracy :math:`\min\left(K,P\right)` if the quadrature can be evaluated with appropriate accuracy. Order reductions may occur if the interpolating polynomial in the quadrature suffers from Runge's phenomenon. As we discuss below, uniformly spaced nodes have some computational advantage but is therefore also associated with some risk. Safer choices include Lobatto nodes or Chebyshev nodes (with inclusion of endpoints) to minimize the risk of order reductions. Implications on the choice of quadrature nodes can be found in :cite:`Layton2005`. 

For explicit advection, the deferred correction procedure integrates the correction equation sequentially and therefore does not allow each substep :math:`\Delta t_m` to exceed the CFL-limited time step :math:`\Delta t_{\textrm{cfl}}`, i.e. :math:`\Delta t_m < \Delta t_{\textrm{cfl}} \forall m`. Since we have :math:`\Delta t = \sum_m\Delta t_m`, uniform nodes maximize :math:`\Delta t` subject to the CFL constraint. For example, an :math:`\verb|SISDC|_P^K` method with uniformly spaced nodes has a maximum possible time step :math:`\Delta t < (P-1)\Delta t_{\textrm{cfl}}`. For the same number of function evaluations, the allowed time step with for Lobatto or Chebyshev nodes is smaller. For :math:`P\leq 3`, the uniform nodes, Lobatto nodes, and Chebyshev nodes coincide. Larger time steps are possible with uniform nodes for :math:`P>3`, which has some computational consequence. The table below summarizes the largest possible time steps for the :math:`\verb|SISDC|_P^K` method with the various quadratures. Finally, note that :math:`\Delta t_m < \Delta t_{\textrm{cfl}}` does not guarantee stability since further restrictions may required for stability of the reaction terms.

==========  =================================== ====================================   ================================
 :math:`P`   Lobatto                             Chebyshev                             Uniform
==========  =================================== ====================================   ================================
2           :math:`\Delta t_{\textrm{cfl}}`      :math:`\Delta t_{\textrm{cfl}}`       :math:`\Delta t_{\textrm{cfl}}`
3           :math:`2\Delta t_{\textrm{cfl}}`     :math:`2\Delta t_{\textrm{cfl}}`      :math:`2\Delta t_{\textrm{cfl}}`
4           :math:`2.26\Delta t_{\textrm{cfl}}`  :math:`1.73\Delta t_{\textrm{cfl}}`   :math:`3\Delta t_{\textrm{cfl}}`
5           :math:`3.05\Delta t_{\textrm{cfl}}`  :math:`2.82\Delta t_{\textrm{cfl}}`   :math:`4\Delta t_{\textrm{cfl}}`
6           :math:`3.50\Delta t_{\textrm{cfl}}`  :math:`3.29\Delta t_{\textrm{cfl}}`   :math:`5\Delta t_{\textrm{cfl}}`
7           :math:`4.26\Delta t_{\textrm{cfl}}`  :math:`4.36\Delta t_{\textrm{cfl}}`   :math:`6\Delta t_{\textrm{cfl}}`
==========  =================================== ====================================   ================================

For the predictor step, it is necessary to evaluate :math:`\mathcal{F}_{\textrm{AR}}\left(\phi_m^{k+1}\right)` and thus update the Poisson and radiative transfer equations at each node. In addition, it is necessary to solve the diffusion equation at every node except :math:`m=0` for every diffusive species, which may also require auxiliary updates of the electric field. The corrector step contains extra floating point operator due to the extra terms :math:`\mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^k\right)` and :math:`\mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^k\right)` and the quadrature :math:`I_m^{m+1}`. The computational cost of adding in these terms is small compared to the cost of an Euler update of the advection-reaction equation since one must also computate source terms, drift velocities, and boundary conditions in addition to construction of the hybrid divergence. In short, the computational cost of the predictor and corrector steps are about the same.

Next, we provide some remarks on the extra computational work involved for higher order methods. Broadly speaking, the total amount of floating point operations increases quadratically with the order. Each node requires evaluation of one advection-reaction operator, at least one electric field update, and one radiative transfer update. Likewise, each substep requires one diffusion solve. Thus, :math:`\verb|SISDC|_K^K` requires :math:`K^2` advection-reaction evaluations, :math:`(K-1)^2` diffusion solves, :math:`(K-1)^2` radiative transfer updates, and at least :math:`K^2` electric field updates. In these estimates we have assumed that the diffusion solve couples semi-implicitly to the electric field, thus each corrector iteration requires one electric field update per node, giving a total cost :math:`K^2`. Strictly speaking, the number of advection-reaction evaluations is slightly less since :math:`\mathcal{F}_{\textrm{AR}}\left(t_0, \phi_0^k\right)` does not require re-evaluation in the corrector, and :math:`\mathcal{F}_{\textrm{AR}}\left(t_p,\phi_p^{K-1}\right)` does not need to be computed for the final iteration since the lagged quadrature is not further needed. Nonetheless, the computational work is quadratically increasing, but this is partially compensated by allowance of larger time steps since the :math:`\verb|SISDC|_K^K` has a stability limit of :math:`(K-1)\Delta t_{\textrm{cfl}}` rather than :math:`\Delta t_{\textrm{cfl}}` for uniformly spaced nodes. For comparison with the predictor :math:`\verb|SISDC|_K^1` which is a first order method, the work done for integration over :math:`(K-1)\Delta t_{\textrm{cfl}}` amounts to :math:`K-1` advection-reaction updates, :math:`K-1` diffusion updates, :math:`K-1` radiative transfer updates, and :math:`K` electric field updates. If we take the electric field updates as a reasonable metric for the computational work, the efficiency of the :math:`K` th order method over the first order method is about :math:`K` for integration over the same time interval, i.e. it increases linearly rather than quadratically. However, this estimate is only valid if we do not take accuracy into account. In practice, the predictor does not provide the same accuracy as the corrector over the same integration interval. A fair comparison of the extra computational work involved would require that the accuracy of the two methods be the same after integration over a time :math:`(K-1)\Delta t_{\textrm{cfl}}`, which will generally require more substeps for the first order method. While we do not further pursue this quantification in this paper, the pertinent point is that the extra computational work involved for tolerance-bound higher order discretizations increases sub-linearly rather than quadratically when compared to lower-order equivalents.

We have implemented the SISDC algorithm in the ``sisdc`` class in :file:`/time_steppers/sisdc`. The following class options are available:

.. literalinclude:: links/sisdc.options

.. _Chap:MISDC:


Stochastic integrators
----------------------

Deterministic CFD integrators are generally not suitable for stochastic ODEs. For example, the recursive nature of Heun's method or spectral deferred corrections hardly make sense for stochastic ODEs.

For fluctuating hydrodynamics we are preparing several temporal integrators. ``godunov`` is generally useful for FHD but we also provide an Euler-Maruyama integrator (``euler_maruyama``) which is first order accurate in time (although with an accurate advective integrator). Although ``euler_maruyama`` is functionally similar to ``godunov``, it does not eliminate the dielectric relaxation time and is therefore less stable for some simulation cases. 

.. _Chap:euler_maruyama:

euler_maruyama
______________

:ref:`Chap:euler_maruyama` implements the Euler-Maruyama method. This method is based on an Euler method with explicit or implicit diffusion.


  
.. bibliography:: references.bib
