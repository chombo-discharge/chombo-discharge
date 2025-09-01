.. _Chap:MeshODESolver:

Mesh ODE solver
===============

The ``MeshODESolver<N>`` implements a solver for

.. math::

   \frac{\partial \vec{\phi}}{\partial t} = \vec{S},

where :math:`\vec{\phi}` represents :math:`N` unknowns on the mesh, and :math:`\vec{S}` is the corresponding source term.
The class is templated as

.. literalinclude:: ../../../../Source/MeshODESolver/CD_MeshODESolver.H
   :language: c++
   :lines: 22-27
   :dedent: 0

where ``N`` indicates the number of variables stored on the mesh.

``MeshODESolver<N>`` is designed to store ``N`` variables in each grid cell, without any cell-to-cell coupling. 
To instantiate the solver, use the full constructor with reference to :ref:`Chap:AmrMesh`:

.. literalinclude:: ../../../../Source/MeshODESolver/CD_MeshODESolver.H
   :language: c++
   :lines: 40-45
   :dedent: 2

.. tip::
   
   Source code for the ``MeshODESolver<N>`` resides in :file:`Source/MeshODESolver`.
   See `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classMeshODESolver.html>`_ for the full C++ API.

Setting :math:`\vec{\phi}`
--------------------------

Mesh-based
__________

To set :math:`\vec{\phi}`, one can fetch the :math:`N` mesh components from

.. literalinclude:: ../../../../Source/MeshODESolver/CD_MeshODESolver.H
   :language: c++
   :lines: 136-141
   :dedent: 2

This will return the data holder that holds the cell centered data :math:`\vec{\phi}` which the user can iterate through.
See :ref:`Chap:MeshData` for examples.

Analytic function
_________________

One can set :math:`\vec{\phi}(\mathbf{x}) = \vec{f}\left(\mathbf{x}\right)` through the following functions:

.. literalinclude:: ../../../../Source/MeshODESolver/CD_MeshODESolver.H
   :language: c++
   :lines: 98-111
   :dedent: 2

These differ only in that the first function sets a specific component, whereas the second version sets all :math:`N` components. 

Setting :math:`\vec{S}`
-----------------------

Mesh-based
__________

For a general method of setting the source term one can fetch :math:`\vec{S}` through

.. literalinclude:: ../../../../Source/MeshODESolver/CD_MeshODESolver.H
   :language: c++
   :lines: 148-152
   :dedent: 2

This returns a reference to :math:`\vec{S}` which the user can iterate through and set the value in each cell.
See :ref:`Chap:MeshData` for explicit examples.

Spatially dependent
___________________

The source term can also be set on a component-by-component basis using

.. literalinclude:: ../../../../Source/MeshODESolver/CD_MeshODESolver.H
   :language: c++
   :lines: 113-119
   :dedent: 2

As a function of :math:`\vec{\phi}`
___________________________________

In order to compute the source term :math:`\vec{S}` as a function of :math:`\vec{\phi}`, ``MeshOdeSolver<N>`` has a function

.. literalinclude:: ../../../../Source/MeshODESolver/CD_MeshODESolver.H
   :language: c++
   :lines: 30-33, 120-126
   :dedent: 2

which computes the source term :math:`\vec{S}` as a function

.. math::

   \vec{S} = \vec{f}\left(\vec{\phi},t\right).

An example which sets :math:`\vec{S} = \vec{\phi}` is given below

.. code-block:: c++

   auto f = [](const std::array<Real, N>& phi, const Real t) -> std::array<Real, N> {
      const std::array<Real, N> S = phi;

      return S;
   };

   solver.computeRHS(f);

Regridding
----------

When regridding the ``MeshODESolver<N>``, one must first ensure that the mesh data on the old mesh is stored before calling the regrid function:

.. literalinclude:: ../../../../Source/MeshODESolver/CD_MeshODESolver.H
   :language: c++
   :lines: 160-167
   :dedent: 2

This must be done *before* :ref:`Chap:AmrMesh` creates the new grids.
This will store :math:`\vec{\phi}` on the old mesh.
After :ref:`Chap:AmrMesh` has generated the new grids, :math:`\vec{\phi}` can be interpolated onto the new grids by calling

.. literalinclude:: ../../../../Source/MeshODESolver/CD_MeshODESolver.H
   :language: c++
   :lines: 169-177
   :dedent: 2

Users can also choose to turn on/off slope limiters when putting the solution on the new mesh.

.. important::
   
   The source term :math:`\vec{S}` is also allocated on the new mesh, but is not interpolated onto the new grids.
   It must therefore be set by the user after calling the regrid function.

Input options
-------------

Several input options are available for configuring the run-time configuration of ``MeshODESolver<N>``, which are listed below

.. literalinclude:: ../../../../Source/MeshODESolver/CD_MeshODESolver.options
   :caption: Input options for the ``MeshODESolver<N>`` class.
	     All options are run-time configurable.   

I/O
---

The user can add :math:`\vec{\phi}` and :math:`\vec{S}` to output files by specifying these in the input script.
These variables are named

.. code-block:: text

   MeshODESolver.plt_vars = phi rhs

Only ``phi`` and ``rhs`` are recognized as valid arguments.
If choosing to omit output variables for the solver, one can put e.g. ``MeshODESolver.plt_vars = -1``. 

.. note::

   ``MeshODESolver<N>`` checkpoint files only contain :math:`\vec{\phi}`. 
   
