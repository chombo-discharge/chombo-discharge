chombo-discharge
----------------

This is chombo-discharge, a multiphysics code which uses Chombo for plasma
simulations with adaptive mesh refinement (AMR) on embedded boundary grids. 

A modified version of Chombo is distributed together with this code.
chombo-discharge only uses Chombo; it is not affiliated nor endorsed by LBNL.

Documentation
_____________
Documentation is available at https://chombo-discharge.github.io

Getting started
_____________

To clone chombo-discharge, set the environment variable ``$DISCHARGE_HOME`` as follows

export DISCHARGE_HOME=<Location for chombo-discharge>

#. ``Chombo`` is included as a submodule in ``chombo-discharge``.
   To clone ``chombo-discharge`` and the dependency ``Chombo``, do

   .. code-block:: bash
		   
      git clone --recursive git@github.com:chombo-discharge/chombo-discharge.git


Contributing
_____________
We welcome feedback, bug reports, or code contributions. Use the github issue tracker and pull request system for code contributions
See code documentation for coding style and review system. 

License
_______

See LICENSE and Copyright.txt for redistribution rights. 
