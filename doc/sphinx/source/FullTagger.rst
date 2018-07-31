.. _Chap:full_tagger:

full_tagger
_______________

:ref:`Chap:full_tagger` is essentially a more general interface than :ref:`Chap:ultralw_tagger` as it makes available more data for tracer evaluation. The only difference between these two classes is that :ref:`Chap:full_tagger` has a tracer computation routine

.. code-block:: c++

     virtual Vector<Real> tracer(const RealVect&         a_pos,
    			         const Real&             a_time,
				 const Real&             a_dx,
				 const RealVect&         a_E,
				 const Real&             a_min_E,
				 const Real&             a_max_E,
				 const RealVect&         a_grad_E,
				 const Real&             a_min_grad_E,
				 const Real&             a_max_grad_E,
				 const Real&             a_rho,
				 const Real&             a_min_rho,
				 const Real&             a_max_rho,
				 const RealVect&         a_grad_rho,
				 const Real&             a_min_grad_rho,
				 const Real&             a_max_grad_rho,
				 const Vector<Real>&     a_ion_densities,
				 const Vector<Real>&     a_min_ion_densities,
				 const Vector<Real>&     a_max_ion_densities,
				 const Vector<RealVect>& a_ion_gradients,
				 const Vector<Real>&     a_min_ion_gradients,
				 const Vector<Real>&     a_max_ion_gradients,
				 const Vector<Real>&     a_photon_densities,
				 const Vector<Real>&     a_min_photon_densities,
				 const Vector<Real>&     a_max_photon_densities) = 0;

Otherwise, this class is essentially equivalent to :ref:`Chap:ultralw_tagger`. Of course, :ref:`Chap:ultralw_tagger` is a subset of :ref:`Chap:full_tagger`, but it's implementation is necessary since it uses less memory.

     
