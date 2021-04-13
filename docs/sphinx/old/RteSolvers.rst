.. _Chap:RteSolvers:

RTE Solvers
-----------

Currently, the radiative transport equation (RTE) solvers are handled in the diffusion approximation. That is, instead of solving the full RTE, we solve

.. math::

   \kappa \Psi - \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi\right) = \frac{\eta}{c},

where :math:`\kappa` is the absorption coefficient and :math:`\eta` is the source term. The interface to these solvers are discussed in the :ref:`Chap:plasma_kinetics` chapter. There are (slightly) more interfaces in place in order to handle more complex solvers in the future.


Currently, there are only dummy options available for selecting solvers:

.. literalinclude:: links/rte_layout.options

The diffusion RTE is also solved by using multigrid, with the implementation class being :doxy:`eddington_sp1 <classeddington__sp1>`. As with :ref:`Chap:PoissonSolvers`, the user may tune the multigrid parameters

.. literalinclude:: links/eddington_sp1.options

The explanation to the various parameters is found in :ref:`Chap:PoissonSolvers`. 
