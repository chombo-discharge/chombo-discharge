.. _Chap:MeshODESolver:

Mesh ODE solver
===============

The ``MeshODESolver`` implements a solver for

.. math::

   \frac{\partial \vec{\phi}}{\partial t} = \vec{S},

where :math:`\vec{\phi}` represents :math:`N` unknowns on the mesh, and :math:`\vec{S}` is the corresponding source term.
The class is templated as

.. code-block:: c++

   template <size_t N = 1>
   class MeshODESolver;

where ``N`` indicates the number of variables stored on the mesh.

To instantiate the solver, use the full constructor with reference to :ref:`Chap:AmrMesh`:

.. code-block:: c++

   template <size_t N>
   MeshODESolver<N>::MeshODESolver(const RefCountedPtr<AmrMesh>& a_amr) noexcept;

If running dual grid simulations, the corresponding :ref:`Chap:Realm` can be set through the public API, see `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classMeshODESolver.html>`_.

.. note::
   
   Source code for the ``MeshODESolver`` resides in :file:`Source/MeshODESolver`.

Setting :math:`\vec{\phi}`
--------------------------

Mesh-based
__________

To set the initial data in a general way, one can fetch the :math:`N` mesh components from

.. code-block::

   template <size_t N>
   EBAMRCellData&
   MeshODESolver<N>::getPhi() noexcept;

This will return the data holder that holds the cell centered data :math:`\vec{\phi}` which the user can iterate through.

Analytic function
_________________

One can set :math:`\vec{\phi}` by using analytic functions :math:`\vec{\phi}(\mathbf{x}) = \vec{f}\left(\mathbf{x}\right)` through

.. code-block:: c++

   using Func1 = const std::function<Real(const RealVect& a_pos)>;

   template <size_t N>
   using Func2 = const std::function<std::array<Real, N>(const RealVect& a_pos)>;

   template <size_t N>
   void
   MeshODESolver<N>::setPhi(const Func1& a_func1, const size_t a_comp) noexcept;

   template <size_t N>
   void
   MeshODESolver<N>::setPhi(const Func2& a_func2) noexcept;

These differ in the sense that ``Func1``, which is just an alias for a function :math:`f = f\left(\mathbf{x}\right)`, sets the value for a specified component.
The other version that takes ``Func2`` as an argument sets the corresponding values for all components. 

Setting :math:`\vec{S}`
-----------------------

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

The above function will evaluate a function :math:`f(\mathbf{x})` in each cell.

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

Regridding
----------

When regridding the ``MeshODESolver``, one should call

.. code-block:: c++

   template <size_t N>
   void
   MeshODESolver<N>::preRegrid(const int a_lbase, const int a_oldFinestLevel) noexcept;

*before* :ref:`Chap:AmrMesh` creates the new grids.
This will store :math:`\vec{\phi}` on the old mesh.
After :ref:`Chap:AmrMesh` has generated the new grids, :math:`\vec{\phi}` can be interpolated onto the new grids by calling

.. code-block:: c++

   template <size_t N>
   MeshODESolver<N>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) noexcept;

Users can also choose to turn on/off slope limiters when putting the solution on the new mesh, see :ref:`Chap:Regridding`.
The source term :math:`\vec{S}` is also allocated on the new mesh, but is not interpolated onto the new grids.

I/O
---

The user can add :math:`\vec{\phi}` and :math:`\vec{S}` to output files by specifying these in the input script.
These variables are named

.. code-block:: txt

   MeshODESolver.plt_vars = phi rhs

Only ``phi`` and ``rhs`` are recognized as valid arguments.
If choosing to omit output variables for the solver, one can put e.g. ``MeshODESolver.plt_vars = -1``. 

.. note::

   ``MeshODESolver`` checkpoint files only contain :math:`\vec{\phi}`. 
