.. _Chap:Compiling:

Compiling PlasmaC
-----------------

Currently, the entire PlasmaC framework is built into your simulation cases. While this is something that we are working on improving, this means that there is no separate build for the PlasmaC source code and your application files.

Once an application is set up, compiling is done by

.. code-block:: bash

   make -s -j 16 DIM=2 <application_name>

Compiling must be perform from the folder which houses your makefile. 
