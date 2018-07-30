.. _Chap:amr_mesh:

amr_mesh
--------

:ref:`Chap:amr_mesh` is the class that handles almost all spatial operations in PlasmaC. Internally, :ref:`Chap:amr_mesh` contains a bunch of operators that are useful across classes, such as ghost cell interpolation operators, coarsening operators, and stencils for interpolation and extrapolation near the embedded boundaries. :ref:`Chap:amr_mesh` also contains routines for generation and load-balancing of grids based and also contains simple routines for allocation and deallocation of memory. For details, see the :doxy:`Doxygen API <amr_mesh>`.

:ref:`Chap:amr_mesh` is an integral part of PlasmaC, and users will never have the need to modify it unless they are implementing entirely new solvers. The behavior of :ref:`Chap:amr_mesh` is modified through it's available input parameters, listed below:

.. literalinclude:: links/amr_mesh.options     

In the above input parameters, **max_amr_depth** indicates the maximum AMR depth that we will simulate. **coarsest_domain**, which is the number of cells on the coarsest AMR level, *must* be divisible by **blocking_factor**, which is the smallest possible allowed box that will be generated in the grid generation. Likewise, **max_box_size** indicates the largest possible box.

**ref_rat** indicates the refinement factors between levels; the first entry indicates the refinement between the coarsest AMR level and the next. We only support refinement factors of 2 or 4 (you may use mixed refinement factors). This means that your coarsest spatial resolution will be given by the **coarsest_domain**, and the finest resolution given by **coarsest_domain** multiplied recursively by your refinement ratios. Note that the resolution in each spatial direction *must* be the same.

**fill_ratio**, which must be :math:`>0` and :math:`\leq 1` indicates the grid generation effectiveness. For higher values of the fill ratio, grids are smaller and more compact, but more boxes are generated. For lower values, grids tend to be larger, and boxes tend to be more square. Note that if **max_box_size** is equal to **blocking_factor**, the generated grids are essentially octrees.

To add some flexibility in how refinement levels are handled in different stages of evolution, we've added an option **max_sim_depth** which restricts the grid generation to a level equal to or lower than **max_amr_depth**. This options exists because checkpoint/erestart (see :ref:`Chap:RestartingSimulations`) do *not* allow changing the AMR levels since this implies a changed geometry as well.

Users may also change the stencils for data interpolation across refinement boundaries, as well as stencils for extrapolation and interpolation near the embedded boundaries. Typically, we use linear interpolation for ghost cells but quadratic interpolation is supported as well, by changing **ghost_interp**. The flags **stencil_order**, **stencil_type**, and **stencil_radius** control the stencils used for interpolation and extrapolation near the EB. 

If users use features that imply that the grid might have crossing between coarse-fine interfaces and embedded boundaries, he should enable the **ebcf** flag. This is the case, for example, if one uses :ref:`Chap:geo_coarsener` to remove parts of the EB mesh. 

The options **num_ghost** and **eb_ghost** should not be changed since much of our code requires three ghost cells. 

