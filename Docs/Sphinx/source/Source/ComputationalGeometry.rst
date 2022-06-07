.. _Chap:ComputationalGeometry:

ComputationalGeometry
======================

``ComputationalGeometry`` is the class that implements geometries in ``chombo-discharge``.
In principle, geometries consist of electrodes and dielectrics but there are many problems where the actual nature of the EB is irrelevant (such as fluid flow). 

.. tip::

   `ComputationalGeometry C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classComputationalGeometry.html>`_.
   
``ComputationalGeometry`` is *not* an abstract class.
The default implementation is an empty geometry -- i.e. a geometry without any objects, but several geometries are included in :file:`$DISCHARGE_HOME/Geometries`.
Making a non-empty ``ComputationalGeometry`` class requires that you inherit from ``ComputationalGeometry`` and instantiate the following class members:

.. code-block:: c++

   Real m_eps0;
   Vector<Electrode> m_electrodes;
   Vector<Dielectric> m_dielectrics;

Here, ``m_eps0`` is the relative gas permittivity, ``m_electrodes`` are the electrodes for the geometry and ``m_dielectrics`` are the dielectrics for the geometry.
These are described in detail further down.

When geometries are created, the ``ComputationalGeometry`` class will first create the (approximations to the) signed distance functions that describe two possible material phases (gas and solid).
Here, the *solid* phase is the part of the computational domain inside the dielectrics, while the *gas phase* is the part of the computational domain that is outside both the electrodes and the dielectrics.

.. _Chap:Electrode:

Electrode
---------

The ``Electrode`` class is responsible for describing an electrode and potentially also its boundary condition.
Internally, this class is lightweight and consists only of a tuple that holds a level-set function and an associated boolean value that tells whether or not the level-set function is at a live voltage or not.
The constructor for the electrode class is:

.. code-block:: c++
   
  Electrode(RefCountedPtr<BaseIF> a_baseIF, bool a_live, Real a_fraction = 1.0);

where the ``a_baseIF`` argument is the level-set function and the ``a_live`` argument is used by some solvers in order to determine if the electrode is at a live voltage or not.
If ``a_live`` is set to false, :ref:`Chap:FieldSolver` will fetch this value and determine that the electrode is at ground. 
Otherwise, if ``a_live`` is set to true then :ref:`Chap:FieldSolver` will determine that the electrode is at live voltage, and the ``a_fraction`` argument is an optional argument that allows the user to set the potential to a specified fraction of the live voltage.

.. tip::

   `Electrode C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classElectrode.html>`_ 

.. _Chap:Dielectric:

Dielectric
----------

The ``Dielectric`` class describes a dielectric.
This class is lightweight and consists of a tuple that holds a level-set function and the associated permittivity.
The constructors are

.. code-block:: c++
   
  Dielectric(RefCountedPtr<BaseIF> a_baseIF, Real a_permittivity);
  Dielectric(RefCountedPtr<BaseIF> a_baseIF, Real (*a_permittivity)(const RealVect a_pos);

where the ``a_baseIF`` argument is the level-set function and the second argument sets a the permittivity.
The permittivity can be set to a constant (first constructor) or to a spatially varying value (second constructor). 

.. tip::

   `Dielectric C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classDielectric.html>`_

Retrieving parts
----------------

It is possible to retrieve the implicit functions for the electrodes and dielectrics through the following member functions:

.. code-block:: c++

   const Vector<Dielectric>& getDielectrics() const;
   const Vector<Electrode>& getElectrodes() const;

Retrieving implicit functions
-----------------------------

When generating the geometry we compute the implicit functions for each *phase*, the gas-phase and the dielectric-phase.
If none of the implicit functions for the electrodes/dielectrics overlap, the resulting function will also be a signed distance function.
Note that these functions are the unions/intersections of all electrodes and dielectrics. 

To retrive the implicit function corresponding to a particular phase, use

.. code-block:: c++

   const RefCountedPtr<BaseIF>& getImplicitFunction(const phase::which_phase a_phase) const;

where ``a_phase`` will be ``phase::gas`` or ``phase::solid``. 
