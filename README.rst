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

.. code-block:: text
		
   export DISCHARGE_HOME=<Location for chombo-discharge>

There are two ways of cloning ``chombo-discharge``. 

#. By including ``Chombo`` as a submodule in ``chombo-discharge``.
   To clone ``chombo-discharge`` and the dependency ``Chombo``, do

   .. code-block:: text
		   
      git clone --recursive git@github.com:chombo-discharge/chombo-discharge.git ${DISCHARGE_HOME}

   Next, set the ``Chombo`` environment variable ``$CHOMBO_HOME`` to ``$DISCHARGE_HOME/Chombo-3.3/lib``, i.e.

   .. code-block:: text

      export CHOMBO_HOME=$DISCHARGE_HOME/Chombo-3.3/lib

#. By cloning ``Chombo`` separately.
   First clone ``chombo-discharge`` without submodules by

   .. code-block:: text
		   
      git clone git@github.com:chombo-discharge/chombo-discharge.git ${DISCHARGE_HOME}

   Next, set the ``$CHOMBO_HOME`` environment variable and clone ``Chombo`` there, i.e.

   .. code-block:: text

      export DISCHARGE_HOME=<Location for chombo-discharge>
      git clone git@github.com:chombo-discharge/Chombo-3-3.git ${CHOMBO_HOME}      
		   


Contributing
_____________
We welcome feedback, bug reports, or code contributions. Use the github issue tracker and pull request system for code contributions
See code documentation for coding style and review system. 

License
_______

See LICENSE and Copyright.txt for redistribution rights. 
