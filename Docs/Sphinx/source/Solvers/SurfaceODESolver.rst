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

.. literalinclude:: ../../../../Source/SurfaceODESolver/CD_SurfaceODESolver.H
   :language: c++
   :lines: 22-27
   :dedent: 0

where :math:`N` indicates the number of variables stored in each cut cell.
This solver is analogous to :ref:`Chap:MeshODESolver`, with the exception that variables are only stored on cut-cells. 

Instantiation
-------------

To instantiate the solver, use one of the following constructors:

.. literalinclude:: ../../../../Source/SurfaceODESolver/CD_SurfaceODESolver.H
   :language: c++
   :lines: 30-43
   :dedent: 2

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

.. literalinclude:: ../../../../Source/SurfaceODESolver/CD_SurfaceODESolver.H
   :language: c++
   :lines: 210-215
   :dedent: 2

This returns a reference to the underlying data which is defined on all cut-cells.
The user can then iterate through this data and set the values accordingly, see :ref:`Chap:MeshIteration`.

Constant value
______________

To set the data directly, ``SurfaceODESolver<N>`` defines functions

.. literalinclude:: ../../../../Source/SurfaceODESolver/CD_SurfaceODESolver.H
   :language: c++
   :lines: 188-200
   :dedent: 2

Setting :math:`\vec{F}`
--------------------------

In order to set the right-hand side of the equation, functions exist that are entirely analogous to the function signatures for setting :math:`\vec{\phi}`:

.. literalinclude:: ../../../../Source/SurfaceODESolver/CD_SurfaceODESolver.H
   :language: c++
   :lines: 224-251
   :dedent: 2

Resetting cells
---------------

``SurfaceODESolver<N>`` has functions for setting values in the subset of the cut-cells representing dielectrics or electrodes.
The information about whether or not a cut-cell is on the dielectric or electrode is passed in through :ref:`Chap:ComputationalGeometry`.

The function signatures are

.. literalinclude:: ../../../../Source/SurfaceODESolver/CD_SurfaceODESolver.H
   :language: c++
   :lines: 298-334
   :dedent: 2

Note that one can always call ``SurfaceODESolver<N>::getPhi()`` to iterate over other types of cell subsets and set the values from there.

Regridding
----------

When regridding the ``SurfaceODESolver<N>``, one must first call call

.. literalinclude:: ../../../../Source/SurfaceODESolver/CD_SurfaceODESolver.H
   :language: c++
   :lines: 278-286
   :dedent: 2

This must be done *before* :ref:`Chap:AmrMesh` creates the new grids, and will store :math:`\vec{\phi}` on the old mesh.
After :ref:`Chap:AmrMesh` has generated the new grids, :math:`\vec{\phi}` can be interpolated onto the new grids by calling

.. literalinclude:: ../../../../Source/SurfaceODESolver/CD_SurfaceODESolver.H
   :language: c++
   :lines: 288-296
   :dedent: 2

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
Users can switch between these two methods by specifying the proper configuration option in the configuration file.

Input options
-------------

Several input options are available for configuring the run-time configuration of ``MeshODESolver``, which are listed below

.. literalinclude:: ../../../../Source/SurfaceODESolver/CD_SurfaceODESolver.options
   :caption: Input options for the ``SurfaceODESolver<N>`` class.
	     All options are run-time configurable.   

.. note::

   ``SurfaceODESolver`` checkpoint files only contain :math:`\vec{\phi}`. 
