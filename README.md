chombo-discharge
----------------

This is ``chombo-discharge``, a multiphysics code which uses ``Chombo`` for discharge simulations with adaptive mesh refinement (AMR) on embedded boundary grids.

A modified version of ``Chombo`` is distributed together with this code.
``chombo-discharge`` only uses ``Chombo``; it is not affiliated nor endorsed by LBNL.

Documentation is available as [HTML](https://chombo-discharge.github.io/chombo-discharge/Base/Installation.html) or as a [PDF](https://github.com/chombo-discharge/chombo-discharge/raw/gh-pages/chombo-discharge.pdf).
Click [here](https://chombo-discharge.github.io/chombo-discharge/Base/Installation.html) for installation instructions. 

License
-------

See LICENSE and Copyright.txt for redistribution rights. 

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
