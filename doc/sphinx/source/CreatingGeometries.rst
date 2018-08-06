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

PlasmaC geometries
------------------

In PlasmaC, all geometries are created by providing ``BaseIF`` instances. For example, a dielectric is defined by an internal class ``dielectric`` (see :doxy:`Doxygen API <dielectric>` for details), and electrodes are supported through the ``electrode`` class (see :doxy:`Doxygen API <electrode>` for details).

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
