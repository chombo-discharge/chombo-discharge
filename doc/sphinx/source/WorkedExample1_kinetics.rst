.. _Chap:WorkedExample1_kinetics:

Physics module
--------------

Here, we define the appropriate physics module that is used for the convection and diffusion of the species. The code segment below defines the entire kinetic framework; we will explain the various implementation functions below.

.. literalinclude:: links/worked_example1_kinetics.H
   :language: c++
		

Advected species
________________

First, we declare the advected species ``my_species``as a nested C++ class within ``my_kinetics``. Note that your species class doesn't not *need* to be nested inside your physics module, but this is often preferred since there might, otherwise, be nameclashes with external code (typically, multiple implementations of an ``electron`` will exist). In the species constructor, we specify the name of the advected species in ``m_name``, and whether or not the species is mobile and diffusive in ``m_mobile`` and ``m_diffusive``. 

.. literalinclude:: links/worked_example1_kinetics.H
   :language: c++
   :lines: 78-107

We have also added input parameters called ``my_kinetics.width`` and ``my_kinetics.center``, with the latter being a vector. Note that ``my_kinetics.center`` must be obtained through array functions. Here, we first create an array holding the spatial coordinates, and then transform this into a vector:

.. literalinclude:: links/worked_example1_kinetics.H
   :language: c++
   :lines: 88-92

Finally, the remaining piece of code returns the initial data for the species. Here, it simply returns a Gaussian

.. math::

   \exp\left(-\frac{(\mathbf{x}-\mathbf{x_0})^2}{2L^2}\right)

This is implemented as

.. literalinclude:: links/worked_example1_kinetics.H
   :language: c++
   :lines: 100-104

Core plasma functions
_____________________

Next, in the constructor for the core kinetics, we simply declare that we wish to use one CDR solver through ``m_num_species``, no RTE solvers through ``m_num_photons``, and then instantiate the species. It is through ``m_species``, which must be the same length as ``m_num_species``, that PlasmaC will instantiate solvers. A single CDR solver is instantiated for each entry in this vector.

.. literalinclude:: links/worked_example1_kinetics.H
   :language: c++
   :lines: 7-12

Next, the function ``compute_cdr_diffusion_coefficients`` is the function responsible for computing the diffusion coefficients and expects that the user returns the diffusion coefficients for all CDR solvers in the order given by ``m_species``. Here, the species is non-diffusive so the return value is arbitrary. 

.. literalinclude:: links/worked_example1_kinetics.H
   :language: c++
   :lines: 16-22

Likewise, the next section of code sets the advection velocity :math:`\mathbf{v} = \mathbf{E}`

.. literalinclude:: links/worked_example1_kinetics.H
   :language: c++
   :lines: 23-29
