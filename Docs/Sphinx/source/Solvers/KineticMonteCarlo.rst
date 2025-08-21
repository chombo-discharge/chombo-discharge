.. _Chap:KineticMonteCarlo:

Kinetic Monte Carlo
===================

Kinetic Monte Carlo (KMC) algoritms are composed of various methods for stochastically simulating chemically reacting systems.
While various flavors of KMC are encountered in different fields of science, KMC in the context of ``chombo-discharge`` is primarily associated with chemistry kernels.

Concept
-------

In ``chombo-discharge`` the Kinetic Monte Carlo solver advances a state (or multiple states) represented e.g. as state vectors

.. math::

   \vec{X}(t) = \begin{pmatrix}
   X_1(t) \\
   X_2(t) \\
   X_3(t) \\
   \vdots
   \end{pmatrix}

Each row in :math:`\vec{X}` represents the population of some chemical species, and is indexed by an integer.
Reactions between species are represented stoichiometrically as

.. math::

   X_A + X_B + \ldots \xrightarrow{k} X_C + X_D + \ldots,

where :math:`k` is the reaction rate.
Each such reaction is associated with a state change in :math:`\vec{X}`, e.g.

.. math::

   \vec{X}\xrightarrow{r} \vec{X} + \vec{\nu}_r,

where :math:`r` is the reaction type and :math:`\nu_r` is the state change associated with the firing of *one* reaction of type :math:`r`.
A set of such reactions is called the *reaction network* :math:`\vec{R}`.

Propensities :math:`a_r\left(\vec{X}\right)` are defined such that :math:`a_r\left(\vec{X}\right)\textrm{d}t` is the probability that exactly one reaction of type :math:`r` occurs in the time interval :math:`[t, t+\textrm{d}t]`.
For unimolecular reactions of the type

.. math::

   X_A + X_B + \ldots \xrightarrow{k} \emptyset

with :math:`A \neq B \neq \ldots` the propensity function is :math:`k X_A X_B \ldots`.
For bimolecular of the type

.. math::

   2X_A \xrightarrow{k} \emptyset

the propensity is :math:`k \frac{1}{2} X_A(X_A-1)` because there are :math:`\frac{1}{2}X_A(X_A-1)` unique pairs of molecules of type A.
Propensities for higher-order reactions can then be expanded using the binomial theorem.
For example, for a third-order reaction :math:`3X_A\xrightarrow{k} \emptyset` the propensity function is :math:`k\frac{1}{6}X_A(X_A-1)(X_A-2)`. 

Various algorithms can be used for advancing the state :math:`\vec{X}` for an arbitrary reaction network :math:`\vec{R}`.

#. The :ref:`Chap:KMCSSAAdvance` (SSA).
   The SSA is also known as the Gillespie algorithm :cite:`Gillespie1977`, and is an exact stochastic solution to the above problem.
   However, it becomes inefficient as the number of reactions per unit time grows. 
   
#. :ref:`Chap:KMCtauAdvance`, which is an approximation to the SSA which uses Poisson sampling of the underlying reactions. 

#. Hybrid advance that switches between the SSA and tau-leaping in their respective limits, see :ref:`Chap:KMCHybridAdvance`.
   The hybrid algorithm is taken from :cite:t:`Cao2006`, and switches between tau leaping and the SSA in their respective limits.

.. _Chap:KMCSSAAdvance:

Stochastic simulation algorithm
-------------------------------

For the SSA we compute the time until the next reaction by

.. math::

   T = \frac{1}{\sum_{r\in\vec{R}} a_r}\ln\left(\frac{1}{u_1}\right)

where :math:`A = \sum_{r\in\vec{R}} a_r` and :math:`u_1` is a uniformly distributed random variable between :math:`0` and :math:`1`.
The type of reaction that fires is deterimined from

.. math::

   r_c = \textrm{smallest integer satisfying } \sum_{r^\prime = 1}^{r_c} a_{r^\prime} > u_2A,

where and :math:`u_2` is another uniformly distributed random variable between :math:`0` and :math:`1`.
The state is then advanced as

.. math::

   \vec{X}(t+T) = \vec{X}(t) + \vec{\nu}_{r_c}.


.. _Chap:KMCtauAdvance:

Tau leaping
-----------

With tau-leaping the state is advanced over a time :math:`\Delta t` as

.. math::

   \vec{X}\left(t+\Delta t\right) =  \vec{X}\left(t\right) + \sum_{r\in\vec{R}} \vec{\nu}_r\mathcal{P}\left(a_r\left[\vec{X}\left(t\right)\right]\Delta t\right),

   
where :math:`\mathcal{P}` is a Poisson-distributed random variable.
Note that tau leaping may fail to give a thermodynamically valid state, and should thus be used in combination with step rejection.

Tau-leaping variants
____________________

The following forms of tau-leaping are also supported:

#. Midpoint tau-leaping.
#. Poisson random-corrections tau-leaping.
#. Implicit Euler-type tau-leaping. 

These methods can be used either as standalone methods or together with the hybrid algorithm.

.. warning::

   We do not recommend implicit methods for reactive problems.
   The reason for this is that exponential growth follows the equation

   .. math::

      \partial_t X = k X,

   where :math:`k>0` is a growth rate.
   Application of the implicit Euler rule to this system yields

   .. math::

      X^{n+1} = \frac{X^n}{1 - k\Delta t},

   which has a pole at :math:`\Delta t = k^{-1}`, and which is only non-negative for :math:`k\Delta t < 1`.
   Similarly, the discretization can then lead to a large overshoot.   
   In practice, the time step then has to be limited to :math:`\Delta t < k^{-1}`.

   On the other hand, an explicit Euler update would yield

   .. math::

      X^{n+1} = \left(1+k\Delta t\right) X^n,

   which is stable for any :math:`\Delta t` (provided :math:`k` is positive).

.. _Chap:KMCHybridAdvance:

Hybrid algorithm
----------------

The hybrid algorithm is taken from :cite:t:`Cao2006`.
Assume that we wish to integrate over some time :math:`\Delta t`, which proceeds as follows:

#. Let :math:`\tau = 0` be the simulated time within :math:`\Delta t`. 
#. Partition the reaction set :math:`\vec{R}` into *critical* and *non-critical* reactions.
   The critical reactions are defined as the subset of :math:`\vec{R}` that are within :math:`N_{\textrm{crit}}` firings away from exhausting one of its reactants.
   The non-critical reactions are defined as the remaining subset.

#. Compute time steps until the firing of the next critical reaction, and a time step such that the propensities of the non-critical reactions do not change by more than some relative factor :math:`\epsilon`.
   Let these time steps be given by :math:`\Delta \tau_{\textrm{c}}`\ and :math:`\Delta \tau_{\textrm{nc}}`.

#. Select a reactive substep within :math:`\Delta t` from

   .. math::

      \Delta \tau = \min\left[\Delta t - \tau, \min\left(\Delta \tau_{\textrm{c}}, \Delta \tau_{\textrm{nc}}\right)\right]

#. Resolve reactions as follows:

   a. If :math:`\Delta \tau_{\textrm{c}} < \Delta \tau_{\textrm{nc}}` and :math:`\Delta \tau_{\textrm{c}} < \Delta t - \tau` then one critical reaction fires.
      Determine the reaction type using the SSA algorithm.

      Next, advance the state using tau leaping for the non-critical reaction.

   b. Otherwise: No crical reactions fire.
      Advance the state using tau-leapnig for the non-critical reactions only.
      An exception is made if :math:`A\Delta\tau` is smaller than some specified threshold in which case we switch to SSA advancement (which is more efficient in this limit). 

#. Check if :math:`\vec{X}` is a thermodynamically valid state.

   a. If the state is valid, accept it and let :math:`\tau \rightarrow \tau + \Delta\tau`.

   b. If the state is invalid, reject the advancement.
      Let :math:`\Delta\tau_{\textrm{nc}} \rightarrow \Delta \tau_{\textrm{nc}}/2` and return to step 4).

#. If :math:`\tau < \Delta t`, return to step 2.

The :cite:t:`Cao2006` algorithm requires algorithmic specifications as follows:

* The factor :math:`\epsilon` which determines the non-critical time step.
* The factor :math:`N_{\textrm{crit}}` which determines which reactions are critical or not.
* Factors for determining when and how to switch to the SSA-based algorithm in step 5b. 

.. _Chap:KMCSolver:

Implementation
--------------

In ``chombo-discharge``, the KMC solver is implemented as

.. literalinclude:: ../../../../Source/KineticMonteCarlo/CD_KMCSolver.H
   :language: c++
   :lines: 40-57, 73-78
   :dedent: 0

Here, the template parameters are:

* ``R`` is the type of reaction to advance with.
* ``State`` is the state vector that the KMC and reactions will advance.
* ``T`` is the internal floating point or integer representation.

.. tip::

   The ``KMCSolver`` C++ API is found at `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classKMCSolver.html>`_.

``KMCSolver`` is designed to operate with the possibility of separating the solver from the reaction and state types.
Several template constraints exist on the reaction type ``R`` as well as the state type ``State``.

State
_____

The ``State`` representation *must* have a member function

.. literalinclude:: ../../../../Source/KineticMonteCarlo/CD_KMCSingleState.H
   :language: c++
   :lines: 85-90
   :dedent: 2

This function should return true if the state is a valid one (e.g., no negative populations) and false otherwise. 
The functionality is used when using the hybrid advancement algorithm, see :ref:`Chap:KMCHybridAdvance`.

Reaction(s)
___________

The reaction representation ``R`` *must* have the following member functions:

.. literalinclude:: ../../../../Source/KineticMonteCarlo/CD_KMCSingleStateReaction.H
   :language: c++
   :lines: 73-108
   :dedent: 2

These template requirements exist so that users can define their states independent of their reactions.
Likewise, reactions can be defined to operate flexibly on state, and the ``KMCSolver`` can be defined without deep restrictions on the states and reactions that are used. 

Defining states
_______________

State representations ``State`` can be defined quite simply (e.g. just a list of indices).
In the absolute simplest case a state can be defined by maintaining a list of populations like below:

.. code-block:: c++

   class MyState {
   public:
      MyState(const size_t numSpecies) {
         m_populations.resize(numSpecies);
      }

      bool isValidState() const {
         return true;
      }
      
      std::vector<long long> m_populations;
   };

More advanced examples can distinguish between different *modes* of populations, e.g. between species that can only appear on the left/right hand side of the reactions.
See :ref:`Chap:KMCDualState` for such an example.

Defining reactions
__________________

See :ref:`Chap:KMCSolver` for template requirements on state-advancing reactions.
Using ``MyState`` above as an example, a minimal reaction that can advance :math:`A\rightarrow B` with a rate of :math:`k=1` is

.. code-block:: c++

   class MyStateReaction {
   public:

      // List of reactants and products
      MyStateReaction(const size_t a_A, const size_t a_B) {
         m_A = a_A;
         m_B = a_B;	 
      }

      // Compute propensity
      Real propensity(const State& a_state) {
         return a_state[m_A];
      }

      // Never consider these reactions to be "critical"
      long long computeCriticalNumberOfReactions(const Mystate& a_state) {
         return std::numeric_limits<long long>::max();
      }

      // Get a vector/list/deque etc. of the reactant's. <some_container> can be e.g. std::vector<size_t> 
      std::list<size_t> R::getReactants() const {
         return std::list<size_t>{m_A};
      }      

      // Get population
      long long population(const size_t& a_reactant, const MyState& a_state) {
         return a_state.m_populations[a_reactant];
      }

      // Advance state with reaction A -> B
      void advanceState(const MyState& s, const long long& numReactions) const {
         s.populations[m_A] -= numReactions;
         s.populations[m_B] += numReactions;
      }

   protected:
      size_t m_A;
      size_t m_B;	 
   };

Advancement routines
____________________

Many advancement routines for the ``KMCSolver`` are defined internally.
The most general one that uses the hybrid advance is

.. literalinclude:: ../../../../Source/KineticMonteCarlo/CD_KMCSolver.H
   :language: c++
   :lines: 382-392
   :dedent: 2

When using the hybrid algorithm, the user should set the hybrid solver parameters through the function

.. literalinclude:: ../../../../Source/KineticMonteCarlo/CD_KMCSolver.H
   :language: c++
   :lines: 104-119
   :dedent: 2

State and reaction examples
---------------------------

``chombo-discharge`` maintains some states and reaction methods that can be useful when solving problems with ``KMCSolver``.
The following two implementations are currently in use:

#. ``KMCSingleState``, see the `KMCSingleState C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classKMCSingleState.html>`_ and the `KMCSingleStateReaction C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classKMCSingleStateReaction.html>`_.
#. ``KMCDualState``, see the `KMCDualState C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classKMCDualState.html>`_ and the `KMCDualStateReaction C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classKMCDualStateReaction.html>`_.

Verification
------------

Verification tests for ``KMCSolver`` are given in

* :file:`$DISCHARGE_HOME/Exec/Convergence/KineticMonteCarlo/C1`
* :file:`$DISCHARGE_HOME/Exec/Convergence/KineticMonteCarlo/C2`  

C1: Avalanche model
___________________

An electron avalanche model is given in :file:`$DISCHARGE_HOME/Exec/Convergence/KineticMonteCarlo/C1`.
The problem solves for a reaction network

.. math::

   X + \emptyset &\xrightarrow{k_i} X + X + \emptyset \\
   X + \emptyset &\xrightarrow{k_a} \emptyset

In the limit :math:`X\gg 1` the exact solution is

.. math::

   X(t) \approx X(0)\exp\left[(k_i-k_a)t\right].

Figure :numref:`Fig:KineticMonteCarloC1` shows the Kinetic Monte Carlo solution for :math:`k_i = 2k_a = 2` and :math:`X(0) = 10`.

.. _Fig:KineticMonteCarloC1:
.. figure:: /_static/figures/KineticMonteCarloC1.png
   :width: 50%
   :align: center

   Comparison of Kinetic Monte Carlo solution with reaction rate equation for an avalanche-like problem.


C2: Schlögl model
_________________

Solution the Schlögl model are given in :file:`$DISCHARGE_HOME/Exec/Convergence/KineticMonteCarlo/C2`.
For the Schlögl model we solve for a single population :math:`X` with the reactions

.. math::

   B_1 + 2X &\xrightarrow{c_1} 3X, \\
   3X  &\xrightarrow{c_2} B1 + 2X, \\
   B2  &\xrightarrow{c_3} X, \\
   X  &\xrightarrow{c_4} B2.   

The states :math:`B_1` and :math:`B_2` are buffered states with populations that do not change during the reactions. 
Figure :numref:`Fig:KineticMonteCarloC1` shows the Kinetic Monte Carlo solutions for rates

.. math::

   c_1 &= 3\times 10^{-7}, \\
   c_2 &= 10^{-4}, \\
   c_3 &= 10^{-3}, \\
   c_4 &= 3.5

and :math:`B_1 = 10^5`, :math:`B_2 = 2\times 10^5`.
The initial state is :math:`X(0) = 250`.

.. _Fig:KineticMonteCarloC2:
.. figure:: /_static/figures/KineticMonteCarloC2.png
   :width: 50%
   :align: center

   Convergence to bi-stable states for the Schlögl model.

