.. _Chap:ultralw_tagger:

ultralw_tagger
_______________

The first interface is the :doxy:`ultralw_tagger <ultralw__tagger>` class for which the user must implement the following functions

.. code-block:: c++

		  /*!
		     @brief Compute tracer field
		  */
		  virtual Vector<Real> tracer(const RealVect&         a_pos,
			                      const Real&             a_time,
					      const Real&             a_dx,
					      const RealVect&         a_E,
					      const Real&             a_min_E,
					      const Real&             a_max_E,
					      const RealVect&         a_grad_E,
					      const Real&             a_min_grad_E,
					      const Real&             a_max_grad_E) = 0;
					      

		  /*!
		     @brief Coarsen a cell based on a tracer field
		  */
		  virtual bool coarsen_cell(const RealVect&         a_pos,
		                            const Real&             a_time,
					    const Real&             a_dx,
					    const int&              a_lvl,
					    const Vector<Real>&     a_tracer,
					    const Vector<RealVect>& a_grad_tracer) = 0;
		  
		  /*!
		     @brief Refine a cell based on a tracer field
		  */
		  virtual bool refine_cell(const RealVect&         a_pos,
		                           const Real&             a_time,
					   const Real&             a_dx,
					   const int&              a_lvl,
					   const Vector<Real>&     a_tracer,
					   const Vector<RealVect>& a_grad_tracer) = 0;

   
The :doxy:`ultralw_tagger <ultralw__tagger>` does all the cell refinement and coarsening operators based on tracer fields (see :ref:`Chap:cell_tagger` for details) which are based *only* on the electric field. An example implentation of this class is :

.. code-block:: c++

	      class my_tagger : public ultralw_tagger {
   
		my_tagger(){
		   m_num_tracers = 1;
		}

		~my_tagger(){}

		Vector<Real> tracer(const RealVect&         a_pos,
		                    const Real&             a_time,
				    const Real&             a_dx,
				    const RealVect&         a_E,
				    const Real&             a_min_E,
				    const Real&             a_max_E,
				    const RealVect&         a_grad_E,
				    const Real&             a_min_grad_E,
				    const Real&             a_max_grad_E){
		   Vector<Real> tracers(1);
		   tracers[0] = a_E.vectorLength();
		}

		bool coarsen_cell(const RealVect&         a_pos,
		                  const Real&             a_time,
				  const Real&             a_dx,
				  const int&              a_lvl,
				  const Vector<Real>&     a_tracer,
				  const Vector<RealVect>& a_grad_tracer){
		   return (a_tracer[0] < 1.0) ? true : false;
		}

		bool refine_cell(const RealVect&         a_pos,
		                 const Real&             a_time,
				 const Real&             a_dx,
				 const int&              a_lvl,
				 const Vector<Real>&     a_tracer,
				 const Vector<RealVect>& a_grad_tracer) {
   		   return (a_tracer[0] > 1.0) ? true : false;
		}
	      };

This fully functional refinement class defined a tracer field :math:`T = |\mathbf{E}|` and refines if :math:`T > 1.0` and coarsens if :math:`T < 0`. One may, of course, define an arbitrary amount of tracer fields. 
