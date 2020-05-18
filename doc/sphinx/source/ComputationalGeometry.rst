.. _Chap:computational_geometry:

computational_geometry
======================

:ref:`Chap:computational_geometry` is the class that implements that geometry.
In `PlasmaC`, we use level-set functions for description of surfaces. :ref:`Chap:computational_geometry` is not an abstract class;
if you pass in an instance of :ref:`Chap:computational_geometry` (rather than a casted instance), you will get a regular geometry without any boundaries.
A new :ref:`Chap:computational_geometry` class requires that you set the following class members:

.. code-block:: c++

   Real m_eps0;
   Vector<electrode> m_electrodes;
   Vector<dielectric> m_dielectrics;

Here, ``m_eps0`` is the gas permittivity, ``m_electrodes`` are the electrodes for the geometry and ``m_dielectrics`` are the dielectrics for the geometry. 

.. _Chap:electrode:

electrode
---------

The :ref:`Chap:electrode` class is responsible for describing an electrode and its boundary conditions. Internally, this class is lightweight and consists only of a tuple that holds a level-set function and an associated boolean value that tells whether or not the level-set function has a live potential or not. The constructor for the electrode class is:

.. code-block:: c++
   
  electrode(RefCountedPtr<BaseIF> a_baseif, bool a_live, Real a_fraction = 1.0);

where the first argument is the level-set function and the second argument is responsible for setting the potential. The third argument is an optional argument that allows the user to set the potential to a specified fraction of the applied potential.

.. _Chap:dielectric:

dielectric
----------

The :ref:`Chap:dielectric` class describes a electrode, this class is lightweight and consists only of a tuple that holds a level-set function and the permittivity.
The constructors are

.. code-block:: c++
   
  dielectric(RefCountedPtr<BaseIF> a_baseif, Real a_permittivity);
  dielectric(RefCountedPtr<BaseIF> a_baseif, Real (*a_permittivity)(const RealVect a_pos);

where the first argument is the level-set function and the second argument sets a constant permittivity (first constructor) or a variable permittivity (second constructor).

How ``computational_geometry`` works
------------------------------------

When geometries are created, the ``computational_geometry`` class will first create the level-set functions that describe two possible material phases (gas and solid) and then compute the cut cell moments.

It is possible to retrieve the level-set functions for each phase, as well as the the electrodes and dielectrics.
This functionality is provided by:

.. code-block:: c++

   const Vector<dielectric>& get_dielectrics() const;
   const Vector<electrode>& get_electrodes() const;

   const RefCountedPtr<BaseIF>& get_gas_if() const;
   const RefCountedPtr<BaseIF>& get_sol_if() const;


