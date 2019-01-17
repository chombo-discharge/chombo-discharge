.. _Chap:NewSimulations:

Setting up simulations
======================

In this section we consider the setup of completely new mini-applications, which is useful for users that want to write their own geometries or plasma kinetics. All exercises use the Python interface, but includes some coding for implementation of physics and geometries. 

Below, :ref:`Chap:WorkedExample1` considers advective motion of a single species. We show how we implement the :ref:`Chap:plasma_kinetics` class in order to instantiate solver and couple it to the boundary conditions. Furthermore, we show how to generate a geometry consisting of two dielectric slabs.

:ref:`Chap:WorkedExample2` uses more complex physics; in this exercise we use an existing implementation of :ref:`Chap:plasma_kinetics` for streamer discharges in air, but we implement a new geometry.

.. toctree::
   :maxdepth: 3

   CreatingGeometries
   WorkedExample1
   WorkedExample2

   
