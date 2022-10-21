.. _Chap:Regridding:

Regridding
==========

Most ``chombo-discharge`` simulations will benefit from using AMR where grids change in time.
When new grids are generated, the following situations occur:

#. Grids were removed, and data underneath the old grids must be replaced by coarsened data from the fine grids.
#. Grids were added, and data on those grids must be interpolated from the coarse data underneath them.

When grids are removed/added, computational particles must be also be associated with the new grids.

Coarsening data
---------------

When removing grids, we simply keep the old grid data.
We thus assume that, prior to regridd, the user coarsened the data appropriately, see :ref:`Chap:Coarsening`.


Interpolation
-------------

``chombo-discharge`` currently supports interpolation onto the new grids using:

#. Constant interpolation, where the fine-grid data is initialized using the value in the coarse grid cells. 
#. With minmod slope-limiters.

The API for interpolating onto the new grids is given in :ref:`Chap:AmrMesh`.
