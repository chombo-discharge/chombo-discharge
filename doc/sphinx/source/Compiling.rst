.. _Chap:Compiling:

Compiling PlasmaC
-----------------

Once an application has been set up, compiling is done by

.. code-block:: bash

   make -s -j 16 DIM=2 <application_name>

Compiling must be performed from the folder which houses your makefile. 

Currently, all of PlasmaC is compiled into your mini-applications. While this is something that we are working on improving, this means that there is no separate build for the PlasmaC source code and your application files.
