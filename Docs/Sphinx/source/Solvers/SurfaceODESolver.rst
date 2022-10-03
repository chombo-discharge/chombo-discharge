.. _Chap:SurfaceODESolver:

Surface ODE solver
==================

``chombo-discharge`` provides a simple solver for ODE equations

.. math::

   \frac{\partial \vec{\phi}}{\partial t} = \vec{F},

where :math:`\vec{\phi}` represents :math:`N` unknowns on the EB.
Note that the underlying data type for :math:`\vec{\phi}` and :math:`\vec{F}` is ``EBAMRIVData``, see :ref:`Chap:MeshData`.
Such a solver is useful, for example, as a surface charge solver where :math:`\phi` is the surface charge density and :math:`F` is the charge flux onto the EB.

The surface charge solver is implemented as

.. code-block:: c++

   template <int N = 1>
   class SurfaceODESolver;

where :math:`N` indicates the number of variables stored in each cut cell.

Instantiation
-------------

To instantiate the solver, use the default constructor

.. code-block:: c++

   template <int N>
   SurfaceODESolver<N>::SurfaceODESolver();

The solver also requires a reference to :ref:`Chap:AmrMesh`, and the computational geometry such that a full instantiation example is

.. code-block:: c++

   SurfaceODESolver<1>* solver = new SurfaceODESolver<1>();

   solver->setAmr(...);
   solver->setComputationalGeometry(...);

Setting :math:`\vec{\phi}`
--------------------------

Mesh-based
__________

To set :math:`\vec{\phi}` on the mesh, one can fetch the underlying data by calling

.. code-block:: c++

   template <int N>
   EBAMRIVData&
   SurfaceODESolver<N>::getPhi() noexcept;

This returns a reference to the underlying data which is defined on all cut-cells.
The user can then iterate through this data and set the values accordingly, see :ref:`Chap:MeshIteration`.

Function-based
______________

To set the data directly, ``SurfaceODESolver`` defines a function

.. code-block:: c++

   template <int N>
   void
   SurfaceODESolver<N>::setPhi(std::function<std::array<Real, N>(const RealVect pos)>& a_func);

where the input argument represents a function :math:`\vec{f} = \vec{f}\left(\mathbf{x}\right)` that returns a value for each component in :math:`\vec{\phi}`.

Setting :math:`\vec{F}`
--------------------------

Mesh-based
__________

To set :math:`\vec{F}` on the mesh, one can fetch the underlying data by calling

.. code-block:: c++

   template <int N>
   EBAMRIVData&
   SurfaceODESolver<N>::getRHS() noexcept;

This returns a reference to the underlying data which is defined on all cut-cells.
The user can then iterate through this data and set the values accordingly, see :ref:`Chap:MeshIteration`.

Function-based
______________

To set the right-hand side directly, ``SurfaceODESolver`` defines a function

.. code-block:: c++

   template <int N>
   void
   SurfaceODESolver<N>::setRHS(std::function<std::array<Real, N>(const RealVect pos)>& a_func);

where the input argument represents a function :math:`\vec{f} = \vec{f}\left(\mathbf{x}\right)` that returns a value for each component in :math:`\vec{F}`.

Resetting cells
---------------

``SurfaceODESolver`` has functions for setting values in the subset of the cut-cells representing dielectrics or electrodes.
The function signatures are

.. code-block::

   template <int N>
   void
   SurfaceODESolver<N>::resetElectrodeCells(const Real a_value);

   template <int N>
   void
   SurfaceODESolver<N>::resetDielectricCells(const Real a_value);

Calling these functions will set the data value in electrode or dielectric cells to ``a_value``.
Note that one can always call ``SurfaceODESolver<N>::getPhi()`` to iterate over other types of cell subsets.

Regridding
----------

When regridding the ``SurfaceODESolver``, one should call

.. code-block:: c++

   template <int N>
   void
   SurfaceODESolver<N>::preRegrid(const int a_lbase, const int a_oldFinestLevel) noexcept;

*before* :ref:`Chap:AmrMesh` creates the new grids.
This will store :math:`\vec{\phi}` on the old mesh.
After :ref:`Chap:AmrMesh` has generated the new grids, :math:`\vec{\phi}` can be interpolated onto the new grids by calling

.. code-block:: c++

   template <int N>
   SurfaceODESolver<N>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) noexcept;

Note that when interpolating to the new grids one can choose to initialize data in the new cells using the value in the underlying coarse cells, i.e.

.. math::

   \vec{\phi}_{\mathbf{i}_{\textrm{fine}}} = \vec{\phi}_{\mathbf{i}_{\textrm{coarse}}}
   
Alternatively one can initialize the fine-grid data such that the area-weighted value of :math:`\vec{\phi}` is conserved, i.e.

.. math::

   \sum_{\mathbf{i}_{\textrm{fine}}}\alpha_{\mathbf{i}_{\textrm{fine}}}\Delta x_{\textrm{fine}}^{D-1}\vec{\phi}_{\mathbf{i}_{\textrm{fine}}} = \alpha_{\mathbf{i}_{\textrm{coar}}}\Delta x_{\textrm{coar}}^{D-1}\vec{\phi}_{\mathbf{i}_{\textrm{coar}}}

which gives

.. math::
   
   \vec{\phi}_{\mathbf{i}_{\textrm{fine}}} = r^{D-1}\frac{\alpha_{\mathbf{i}_{\textrm{coar}}}}{\sum_{\mathbf{i}_{\textrm{fine}}}\alpha_{\mathbf{i}_{\textrm{fine}}}}\vec{\phi}_{\mathbf{i}_{\textrm{coar}}},

where :math:`\mathbf{i}_{\textrm{fine}}` is set of cut-cells that occur when refining the coarse-grid cut-cell :math:`\mathbf{i}_{\textrm{coar}}` and :math:`r` is the refinement factor between the two grid levels. 
In this case :math:`\vec{\phi}` is strictly conserved.
Users can switch between these two methods by specifying the regridding type in the input script:

.. code-block:: text

   SurfaceODESolver.regridding = arithmetic

or

.. code-block:: text

   SurfaceODESolver.regridding = conservative

I/O
---

The user can add :math:`\vec{\phi}` and :math:`\vec{F}` to output files by specifying these in the input script.
These variables are named

.. code-block:: text

   SurfaceODESolver.plt_vars = phi rhs

Only ``phi`` and ``rhs`` are recognized as valid arguments.
If choosing to omit output variables for the solver, one can put e.g. ``SurfaceODESolver.plt_vars = -1``. 

.. note::

   ``SurfaceODESolver`` checkpoint files only contain :math:`\vec{\phi}`. 
