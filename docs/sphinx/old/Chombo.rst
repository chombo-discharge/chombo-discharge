.. _Chap:Chombo:

Chombo coding guide
-------------------

PlasmaC is, in essence, a large `Chombo <https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations>`_ application. Chombo uses dimensionless structures and wappers that the user should familiarize himself with. The most important structures are

* :file:`Real` - a replacement for float or double (depending on your compiler settings)
* :file:`RealVect` - a vector in space.
* :file:`Vector` - a wrapper for :file:`std::vector`.
* :file:`RefCountedPtr<T>` - a pointer class with reference counting and auto-deallocation.

The useage of these classes is straightforward. For example, a :file:`Real` is declared

.. code-block:: c++

		Real foo = 1.0;

:file:`RealVect` is a spatial vector that contains two or three entries in PlasmaC. To use :file:`RealVect`, one may do

.. code-block:: c++

		RealVect foo = RealVect(1.0, 0.0);

in two dimensions and

.. code-block:: c++

		RealVect foo = RealVect(1.0, 0.0, 0.0);

in three dimensions. The dimensionless way of doing this is to use Chombo macros; 

.. code-block:: c++

		RealVect foo = RealVect(D_DECL(1.0, 0.0, 0.0));

where :file:`D_DECL` is macro that returns the first two variables in 2D, and all three variables in 3D.

The :file:`Vector` class is used just as :file:`std::vector`: 		

.. code-block:: c++

   Vector<Real> foo(2);
   foo[0] = 1.0;
   foo[1] = 0.0;
		
The same goes with the smart pointer :file:`RefCountedPtr<T>`:
   
.. code-block:: c++

   RefCountedPtr<Real> ptr = RefCountedPtr<Real> (new Real(0.0));

For the full Chombo API, please see the `Chombo doxygen guide <http://davis.lbl.gov/Manuals/CHOMBO-RELEASE-3.2/classes.html>`_. 
