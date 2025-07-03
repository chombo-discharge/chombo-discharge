.. _Chap:ComputationalGeometry:

ComputationalGeometry
=====================

``ComputationalGeometry`` is the class that implements geometries in ``chombo-discharge``.
In principle, geometries consist of electrodes and dielectrics but there are many problems where the actual nature of the EB is irrelevant (such as fluid flow).
For other problems, such as ones involving electric fields, the classification into electrodes and dielectric are obviously important.

.. tip::

   Here is the `ComputationalGeometry C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classComputationalGeometry.html>`_.
   
``ComputationalGeometry`` is *not* an abstract class.
The default implementation is an empty geometry, i.e., a geometry without any solid objects.
Several pre-defined geometries are included in :file:`$DISCHARGE_HOME/Geometries`.

Making a non-empty ``ComputationalGeometry`` class requires that you inherit from ``ComputationalGeometry`` and instantiate the class members specified in :numref:`ComputationalGeometryMembers`.

.. _ComputationalGeometryMembers:		    
.. literalinclude:: ../../../../Source/Geometry/CD_ComputationalGeometry.H
   :language: c++
   :caption: List of all data member of the ``ComputationalGeometry`` base class. Highlighted members must be instantiated by the user in order to create a new geometry.
   :lines: 158-196
   :emphasize-lines: 4, 24, 29	   
   :dedent: 2

Here, ``m_eps0`` is the *relative gas permittivity* (which is virtually always 1), ``m_electrodes`` are the electrodes for the geometry and ``m_dielectrics`` are the dielectrics for the geometry.
These are described in detail below.

Implicit functions and compound geometries
------------------------------------------

``ComputationalGeometry`` always uses implicit functions for describing the geometry, see :ref:`Chap:GeometryRepresentation`.
These functions are analytic functions that describe the inside and outside regions of a solid object, and usually appear in the form :math:`f: \mathbb{R}^3\rightarrow \mathbb{R}`.
We point out that there *is* support in ``chombo-discharge`` for turning surface meshes into such functions, and arbitrary complex geometries can therefore be generated.
In ``chombo-discharge`` we always use signed distance functions wherever we can. 

When geometries are created, the ``ComputationalGeometry`` class will first create the (approximations to the) compound signed distance functions that describe two possible material phases (gas and solid).
Here, the *solid* phase is the part of the computational domain inside the dielectrics, while the *gas phase* is the part of the computational domain that is outside both the electrodes and the dielectrics.

Because of this, there is a corresponding logic that underpins the partition into the three relevant domains.
Let :math:`A\cup B` be the union of two objects :math:`A` and :math:`B`, and :math:`A \cap B` be the intersection of the same objects.
Furthermore, let :math:`f_1, f_2, \ldots` denote the implicit functions for the electrodes.
The electrode region is then given by


.. math::    
   :label: ElectrodeRegion

   f = f_1 \cup f_2 \cup \ldots.	   



For the dielectric region, the canonical definition in ``chombo-discharge`` is that this region is composed of all regions that are inside a dielectric but outside of a provided electrode (note that the user can override this behavior in his/her implementation of ``ComputationalGeometry``).
Specifically, if :math:`g_1, g_2, \ldots` are the dielectric implicit functions, the dielectric region is given by


.. math::
   :label: DielectricRegion

   g = \left(g_1 \cup g_2 \cup \ldots \right) \cap f^\complement

The region occupied by the gas is then given by :math:`\left(f \cup g\right)^\complement`.

.. _Chap:Electrode:

Electrode
---------

The ``Electrode`` class is responsible for describing an electrode and also its boundary electrostatic boundary condition.
Internally, this class is lightweight and consists only of a tuple that holds an implicit function and an associated boolean value that tells whether or not the level-set function is at a live voltage or not.
The constructor for the electrode class is given in :numref:`List:Electrode`.

.. _List:Electrode:
.. literalinclude:: ../../../../Source/Geometry/CD_Electrode.H
   :caption: Constructor for the the ``Electrode`` object.
   :lines: 40
   :language: c++
   :dedent: 2

In the code block above, ``a_baseIF`` argument is the implicit function for the electrode object, and the ``a_live`` argument is used by some solvers in order to determine if the electrode is at a live voltage or not.
For example, if ``a_live`` is set to false, :ref:`Chap:FieldSolver` will fetch this value and determine that the electrode is at ground. 
Otherwise, if ``a_live`` is set to true then :ref:`Chap:FieldSolver` will determine that the electrode is at live voltage, and the ``a_fraction`` argument is an optional argument that allows the user to set the potential to a specified fraction of the live voltage.
This is also used throughout other physics modules in ``chombo-discharge``.
The user can also set a specified fraction of the live voltage but setting the ``a_voltageFraction`` parameter to a relative fraction of the applied voltage.

.. tip::

   Here is the `Electrode C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classElectrode.html>`_ 

.. _Chap:Dielectric:

Dielectric
----------

The ``Dielectric`` class describes a dielectric similar to how ``Electrode`` describes an electrode object. 
This class is lightweight and consists of a tuple that holds a level-set function and the relative associated permittivity (just like ``Electrode`` holds an implicit function and a relative voltage). 
The constructors for this class are

.. literalinclude:: ../../../../Source/Geometry/CD_Dielectric.H
   :caption: Constructors for the the ``Dielectric`` object.
   :language: c++
   :lines: 36-50
   :dedent: 2

where the ``a_baseIF`` argument is the level-set function and the second argument sets a the permittivity.
In the constructor, the relative permittivity can be set to a constant (first constructor) or to a spatially varying value (second constructor).
Several solvers will then use the permittivities were specified by the user.

.. tip::

   Here is the  `Dielectric C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classDielectric.html>`_

Retrieving parts
----------------

It is possible to retrieve the implicit functions for the electrodes and dielectrics through the following member functions:

.. literalinclude:: ../../../../Source/Geometry/CD_ComputationalGeometry.H
   :caption: Member functions for retrieving the defined electrodes and dielectrics. 
   :language: c++
   :lines: 48-60
   :dedent: 2

Obtaining the implicit functions for each part can be useful when determining which object is closest to some physical location :math:`\mathbf{x}`.
This is frequently used when, e.g., colliding particles with the embedded boundaries and we want to determine which electrode or dielectric the particle collided with.

Retrieving compound implicit functions
--------------------------------------

When generating the geometry we compute the implicit functions for each *phase*, the gas-phase, the dielectric-phase, and the electrodes.
We outlined this process in :eq:`ElectrodeRegion` and :eq:`DielectricRegion`.

To retrive the implicit function corresponding to a particular phase, use

.. code-block:: c++

   const RefCountedPtr<BaseIF>& getImplicitFunction(const phase::which_phase a_phase) const;

where ``a_phase`` will be ``phase::gas`` or ``phase::solid``.
This function will return the implicit function corresponding to all boundaries that contain the gas phase (``phase::gas``) or dielectric phase (``phase::solid``).
There is no corresponding function for the interior of the electrodes as the solutions are uninteresting in these regions.


