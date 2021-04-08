.. _Chap:Directories:

Directories
-----------

The following directories in PlasmaC are worth noting:

* :file:`/src` contains the PlasmaC source code discussed in :ref:`Chap:ImportantClasses`. 
 
  * :file:`/src/amr_mesh` contains :ref:`Chap:amr_mesh` related code
  * :file:`/src/cdr_solver` contains code for the CDR solvers
  * :file:`/src/elliptic` contains operators for elliptic equations (mostly multifluid Poisson stuff)
  * :file:`/src/geometry` contains code related to the geometric interface
  * :file:`/src/global` contains some globally useful code, such as data structures, stencil types and so on.
  * :file:`/src/plasma_solver` contains the plasma framework, i.e. :ref:`Chap:plasma_kinetics`, :ref:`Chap:plasma_engine` and some related code.
  * :file:`/src/poisson_solver` contains the abstract Poisson solver class and it's geometric multigrid implementation.
  * :file:`/src/rte_solver` contains the RTE solvers
  * :file:`/src/sigma_solver` contains the surface charge solver
* :file:`/geometries_prebuilt` contains some predefined geometries.
* :file:`/plasma_models` and its subdirectories contains various implementation of :ref:`Chap:plasma_kinetics`. 
* :file:`/cell_taggers` and its subdirectories contains various implementation of :ref:`Chap:cell_tagger`.
* :file:`/time_steppers` and its subdirectories contains various implementation of :ref:`Chap:time_stepper`.
* :file:`/base_tests` contains some base tests of PlasmaC
* :file:`/doc` contains the documentation of PlasmaC
    
  * :file:`/doc/sphinx` contains the Sphinx documentation
  * :file:`/doc/doxygen` contains some markup used for the :doxy:`Doxygen API <index>`.
  * :file:`/doc/figures` contains some figures used throughout the documentation. 
* :file:`/app_builder` contains the Python interface for setting up mini-applications.


If you want to extend the PlasmaC code, you *may* write your own mini-apps outside of the PlasmaC framework. However, for maximum reuseability you might want to ensure that your changes are available in the future as well. We recommend that you place your geometries, plasma kinetics, and cell taggers in the appropriate directories listed above. This will also ensure that your work can be reached through our :ref:`Chap:PythonInterface`. 
