.. _Chap:MeshODESolver:

Mesh ODE solver
===============

The ``MeshODESolver`` implements a solver for

.. math::

   \frac{\partial \vec{\phi}}{\partial t} = \vec{S},

where :math:`\vec{\phi}` represents :math:`N` unknowns on the mesh, and :math:`\vec{S}` is the corresponding source term.
The class is templated as

.. code-block:: c++

   template<size_t N = 1>
   class MeshODESolver;

where ``N`` indicates the number of variables stored on the mesh. 

Initializing MeshODESolver
--------------------------

Setting initial data
--------------------

Setting source terms
--------------------

General approach
________________

For a general method of setting the source term one can fetch :math:`\vec{S}` through

.. code-block:: c++

   template <size_t N>
   EBAMRCellData&
   MeshODESolver<N>::getRHS() noexcept;

This returns an l-value reference to :math:`\vec{S}` which the user can iterate through and set the value in each cell.

Spatially dependent
___________________

The source term can also be set on a component-by-component basis using

.. code-block:: c++

   using Func = std::function<Real(const RealVect&)>;

   template <size_t N>
   void
   MeshODESolver<N>::setRHS(const Func& a_rhsFunction, const size_t a_comp) noexcept;

The above function will evaluate a function :math:`S_{\text{a_comp}} = f(\mathbf{x})` in each cell.

Analytically coupled
____________________

A third option is to compute the right-hand side directly using a coupling function.
``MeshODESolver`` aliases a function

.. code-block:: c++

   template<size_t N>
   using RHSFunction = std::function<std::array<Real, N>(const std::array<Real, N>& phi, const Real& time)>;

which computes the source term :math:`\vec{S}` as a function

.. math::

   \vec{S} = \vec{f}\left(\vec{\phi},t\right)

To fill the source term using an analytic coupling function like this, one can call

.. code-block:: c++

   template <size_t N>
   void
   MeshODESolver<N>::computeRHS(const RHSFunction& a_rhsFunction) noexcept;		
