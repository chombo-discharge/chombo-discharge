.. _Chap:photon:

photons
-------

:ref:`Chap:photon` is the class that supplies extra information to the RTE solvers. In those solvers, the source term computation is handled by :ref:`Chap:plasma_kinetics`, so the :ref:`Chap:photon` class is very lightweight. The user must implement a single function which specifies the absorption coefficient at a point in space:

.. code-block:: c++

		virtual const Real get_absorption_coeff(const RealVect& a_pos) const = 0;

In addition, the user should provide a name for the RTE solver so that it can be identified in the output files. This is done by setting a ``m_name`` attribute in the :ref:`Chap:photon` class.

The following is a full implementation of the :ref:`Chap:photon` class:

.. code-block:: c++

		class my_photon : public photon {
		  my_photon() {
		     m_name = "my_photon";
		  }

		  ~my_photon(){}

		  const Real get_absorption_coeff(const RealVect& a_pos) const {
		     return 1.0;
		  }
		};

By default, there are no input parameters available for the :ref:`Chap:photon` class, but the user will often want to include these, for example by modifying the absorption coefficient. Note that you are allowed to use a spatially varying absorption coefficient. Please see :ref:`Chap:MiniApplications` for how to pass input parameters into your classes. 
