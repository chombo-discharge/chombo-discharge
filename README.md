chombo-discharge
----------------

This is ``chombo-discharge``, a multiphysics code which uses ``Chombo`` for discharge simulations with adaptive mesh refinement (AMR) on embedded boundary grids.

A modified version of ``Chombo`` is distributed together with this code.
``chombo-discharge`` only uses ``Chombo``; it is not affiliated nor endorsed by LBNL.

<img src="./Docs/Sphinx/source/_static/figures/BranchingAir.gif" width="75%">


Installation
------------

A serial build quickstart is given below. 
For complete installation instructions, see https://chombo-discharge.github.io/chombo-discharge/Base/Installation.html


Documentation
-------------

User documentation is available as [HTML](https://chombo-discharge.github.io/chombo-discharge/) or as a [PDF](https://github.com/chombo-discharge/chombo-discharge/raw/gh-pages/chombo-discharge.pdf).
A doxygen-generated API is [also available](https://chombo-discharge.github.io/chombo-discharge/doxygen/html/index.html).

License
-------

See LICENSE and Copyright.txt for redistribution rights.


Serial build quickstart
-----------------------

For doing a quick clone and test build of ``chombo-discharge`` without HDF5 capabilities, execute the following steps:

1. Clone chombo-discharge, including third-party dependencies

   ```
   git clone --recursive git@github.com:chombo-discharge/chombo-discharge.git
   ```

1. Install the LAPACK, BLAS, and GCC dependencies:

   ```
   sudo apt install csh gfortran g++ libblas-dev liblapack-dev
   ```
   
2. Choose an installation directory and clone ``chombo-discharge`` there:

   ```
   export DISCHARGE_HOME=/home/foo/chombo-discharge		
   export CHOMBO_HOME=$DISCHARGE_HOME/Submodules/Chombo-3.3/lib
		
   git clone --recursive git@github.com:chombo-discharge/chombo-discharge.git ${DISCHARGE_HOME}   
   ```

3. Copy the GNU compiler makefile to it's appropriate location

   ```
   cp $DISCHARGE_HOME/Lib/Local/Make.defs.GNU $CHOMBO_HOME/mk/Make.defs.local
   ```

4. Build ``chombo-discharge``

   ```
   cd $DISCHARGE_HOME
   make -s -j4
   ```

5. Run a simple example program

   ```
   cd $DISCHARGE_HOME/Exec/Examples/AdvectionDiffusion/DiagonalFlowNoEB
   make -s -j4
   ./*.ex example.inputs
   ```

6. Run an advanced example program

   ```
   cd $DISCHARGE_HOME/Exec/Examples/CdrPlasma/StochasticAir
   make -s -j4
   ./*.ex positive2d.inputs
   ```		

Contributing
------------

We welcome feedback, bug reports, or code contributions.

1. Create a branch for the new feature.

   ```
   git checkout main
   git pull
   git checkout -b my_branch
   ```
   
2. Develop the feature.

   ```
   git add .
   git commit -m "my commit message"
   ```

   If relevant, add Sphinx and doxygen documentation.
   
3. Format the source and example codes using ```clang-format```:

   ```
   find Source Physics Geometries Exec \( -name "*.H" -o -name "*.cpp" \) -exec clang-format -i {} +
   ```
   
4. Push the changes to GitHub

   ```
   git push --set-upstream origin my_branch
   ```
   
5. Create a pull request and make sure the GitHub continuous integration tests pass.
