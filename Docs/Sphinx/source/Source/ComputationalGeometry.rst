.. _Chap:ComputationalGeometry:

ComputationalGeometry
======================

``ComputationalGeometry`` is the class that implements geometries in ``chombo-discharge``.
The source for for this class is located in :file:`$DISCHARGE_HOME/Source/Geometry/`.

.. note::

   ``ComputationalGeometry`` is *not* an abstract class.
   The default implementation is an empty geometry -- i.e. a geometry without any objects. 

Making a non-empty ``ComputationalGeometry`` class requires that you inherit from ``ComputationalGeometry`` and set the following class members:

.. code-block:: c++

   Real m_eps0;
   Vector<Electrode> m_electrodes;
   Vector<Dielectric> m_dielectrics;

Here, ``m_eps0`` is the relative gas permittivity, ``m_electrodes`` are the electrodes for the geometry and ``m_dielectrics`` are the dielectrics for the geometry.
These are described in detail further down.

When geometries are created, the ``ComputationalGeometry`` class will first create the (approximations to the) signed distance functions that describe two possible material phases (gas and solid).
Here, the *solid* phase is the part of the computational domain inside the dielectrics, while the *gas phase* is the part of the computational domain that is outside the electrodes and dielectrics.

.. _Chap:Electrode:

Electrode
---------

The ``Electrode`` class is responsible for describing an electrode and its boundary conditions. Internally, this class is lightweight and consists only of a tuple that holds a level-set function and an associated boolean value that tells whether or not the level-set function has a live potential or not. The constructor for the electrode class is:

.. code-block:: c++
   
  Electrode(RefCountedPtr<BaseIF> a_baseif, bool a_live, Real a_fraction = 1.0);

where the first argument is the level-set function and the second argument is responsible for setting the potential. The third argument is an optional argument that allows the user to set the potential to a specified fraction of the applied potential.

.. _Chap:Dielectric:

Dielectric
----------

The ``Dielectric`` class describes a dielectric.
This class is lightweight and consists only of a tuple that holds a level-set function and the associated permittivity.
The constructors are

.. code-block:: c++
   
  Dielectric(RefCountedPtr<BaseIF> a_baseif, Real a_permittivity);
  Dielectric(RefCountedPtr<BaseIF> a_baseif, Real (*a_permittivity)(const RealVect a_pos);

where the first argument is the level-set function and the second argument sets a constant permittivity (first constructor) or a variable permittivity (second constructor).


Retrieving distance functions
-----------------------------

It is possible to retrieve the SDFs for each phase, as well as the the electrodes and dielectrics.
This functionality is provided by:

.. code-block:: c++

   const Vector<Dielectric>& getDielectrics() const;
   const Vector<Electrode>& getElectrodes() const;

   const RefCountedPtr<BaseIF>& getGasImplicitFunction() const;
   const RefCountedPtr<BaseIF>& getSolidImplicitFunction() const;


