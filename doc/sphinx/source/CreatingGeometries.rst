.. _Chap:CreatingGeometries:

Creating geometries
===================

In PlasmaC, support for internal boundary is done by means of an embedded boundary approach. The surface itself is computed by intersecting the mesh with a level-set function, i.e. a function of the form

.. math::

   f\left(\mathbf{x}\right) = 0,

where :math:`f\left(\mathbf{x}\right) > 0` is the gas phase, :math:`f\left(\mathbf{x}\right) = 0` defines the interface, and :math:`f\left(\mathbf{x}\right) < 0` defines the solid phase.


BaseIF
------

In Chombo, level-set functions are represented by a ``BaseIF`` object that the user must implement. There is native support for a wide variety of shapes in Chombo, as well as tools for performing constructive solid geometry (CSG).

The user *must* implement two functions for creating a new geometric shape:

.. code-block:: c++
		
   virtual Real value(const RealVect& a_point) const = 0

and

.. code-block:: c++
		
   virtual BaseIF* newImplicitFunction() const = 0

where the latter method is a factory method that produces a new instance of the derived class. The Chombo library has multiple examples of how ``BaseIF`` is used. If you don't know how to write a ``BaseIF`` derived class, you should check out the Chombo implementation of some of the basic shapes, or some of the basic shapes below.

.. rubric:: Implemented ``BaseIF`` classes

We have various implemented of ``BaseIF`` that the user may or may not find useful:

* ``box_if`` A 2D or 3D box with rounded corners.
* ``cylinder_if`` A cylindrical shape
* ``dcel_if`` A polygonal tessellation - i.e. support for triangulated objects
* ``graded_perlin_sphere_if`` A sphere with surface noise (roughness) on one side
* ``hyperboloid_if`` A hyperboloid
* ``perlin_if`` Implemented of Ken Perlin's landscape noise function
* ``perlin_plane_if`` A plane with surface roughness (implemented with Perlin noise)
* ``perlin_rod_if`` A cylindrical rod with surface roughness at the end cap
* ``perlin_sphere_if`` A sphere with surface roughness
* ``rod_if.H`` A needle rod, essentially a cylinder with a spherical rounded end cap

If you want to use these function, you should refer to the :doxy:`Doxygen API <index>` for the class that you're interested in.

We have extended Chombo with support for polygon tessellations through the ``dcel_if`` class which uses a doubly connected edge list (DCEL) in order to describe a watertight surface mesh. However, the current interface is not particularly user-friendly. If you are interested in using this, please :ref:`Chap:contact`.

.. figure:: figures/car.png
   :width: 480px
   :align: center

   Cut-cell grid generation from surface tessellations

PlasmaC geometries
------------------

In PlasmaC, all geometries are created by providing ``BaseIF`` instances that describe the shape of a solid, and then pass these into the ``dielectric`` and ``electrode`` classes. These two classes are lightweight classes that pass additional information into PlasmaC, for example the permittivity or live voltage of an electrode. 

Dielectrics
___________

A dielectric is created through one of the constructors

.. code-block:: c++

   /*!
      @brief Constructor. 
   */
   dielectric(RefCountedPtr<BaseIF> a_baseif, Real a_permittivity);

   /*!
      @brief Constructor, allows for spatially dependent permittivity
   */
   dielectric(RefCountedPtr<BaseIF> a_baseif, Real (*a_permittivity)(const RealVect a_pos));

To create a ``dielectric``, the user needs to provide the *shape*, i.e. the ``BaseIF`` instance that describes the surface of the dielectric, and a constant or variable permittivity (the second constructor allows for spatially varying permittivities). 

Electrodes
__________

Likewise, electrodes are created through the ``electrode`` class constructor

.. code-block:: c++

   electrode(RefCountedPtr<BaseIF> a_baseif, bool a_live, Real a_fraction = 1.0);

Here, the constructor takes the shape through a ``BaseIF`` instance, and the user can specify if the electrode is at live voltage, and also specify the fraction of the live voltage. 

Example geometry
________________

PlasmaC expects you to instantiate your geometries through a weak constructor. For example, the code segment below defines a geometry consisting only of a single electrode defined as a sphere

.. literalinclude:: example_geometry.cpp
   :language: c++

The above piece of code defined a sphere at live voltage sitting at the origin with radius 1. Creating dielectrics is done in precisely the same way, except that you need to pass in the *permittivity* instead of the live voltage. 
