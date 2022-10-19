.. _Chap:KineticMonteCarlo:

Kinetic Monte Carlo
===================

The 

.. _Chap:KMCSolver:

Implementation
--------------

The Kinetic Monte Carlo solver is implemented as

.. code-block:: c++

   template <typename R, typename State, typename T = long long>
   class KMCSolver
   {
   public:
      using ReactionList = std::vector<std::shared_ptr<const R>>;
      
      inline KMCSolver(const ReactionList& a_reaction) noexcept;
   }

The template parameters are:

* ``R`` is the type of reaction to advance with.
* ``State`` is the state vector that the KMC and reactions will advance.
* ``T`` is the integer representation.

.. tip::

   The ``KMCSolver`` C++ API is found at `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classKMCSolver.html>`_.

State
_____

The ``State`` representation *must* have a member function

.. code-block:: c++

   bool State::isValidState() const;

which determines if the state is thermodynamically valid (e.g. no negative populations).
The functionality is used when using the hybrid advancement algorithm, see :ref:`Chap:KMCHybridAdvance`.

Reaction(s)
___________

The reaction representation ``R`` *must* have the following member functions:

.. code-block:: c++

   // Compute the propensity of the current reaction. 
   Real R::propensity(const State& s) const;

   // Compute the number of reactions before exhausting one of the reactants
   T R::computeCriticalNumberOfReactions(const State& s) const;

   // Compute the number of reactions before exhausting one of the reactants
   void R::advanceState(const State& s, const T& numReactions) const;

   // Get a vector/list/deque etc. of the reactants. <some_container> can be e.g. std::vector<size_t> 
   <some_container> R::getReactants() const;

   // Get the population corresponding to 'reactant' in the input state. If e.g. <some_container> is
   // std::vector<size_t> then <some_type> will be <size_t>
   T R::population(const <some_type> reactant, const State& s) const;

These template requirements exist so that users can define their states independent of their reactions.
Likewise, reactions can be defined to operate flexibly on state, and the ``KMCSolver`` can be defined without deep restrictions on the states and reactions that are used. 

Defining states
---------------

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
------------------

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

State and reaction examples
---------------------------

``chombo-discharge`` 

.. _Chap:KMCSingleState:

Single-state
____________

.. _Chap:KMCDualState:

Dual-state
__________


Advancement algorithms
----------------------

Stochastic simulation algorithm
_______________________________

:math:`\tau` leaping
____________________


.. _Chap:KMCHybridAdvance:

Hybrid algorithm
________________

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

