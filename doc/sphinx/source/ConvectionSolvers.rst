.. _Chap:ConvectionSolvers:

CDR solvers
-----------

In PlasmaC, the convection-diffusion-reaction solvers are represented by type-casted instances of :doxy:`cdr_solver <classcdr__solver>`. This class contains the important function that computes the right hand side :math:`\nabla\cdot(D\nabla n - \mathbf{v}n) + S`, which is the most important function that allows for using the method of lines. It also contains most of the functionality required for writing new discretizations of the convection-diffusion-reaction equations.

Selecting solvers
_________________

In PlasmaC, solvers are selected at run-time by providing appropriate input parameters. Currently, we support the Scharfetter-Gummel scheme and a Godunov method which computes fluxes by finding the slope-limited face states. These solvers are selected by modifying the input parameters:

.. literalinclude:: links/cdr_layout.options

In the above, **which_solver** specifies the solver that user wants. For the Scharfetter-Gummel scheme, there are no further options. For the Godunov methods, the user may also choose the way the hybrid advective divergence is computed in boundary cells by modifying **cdr_gdnv.divF_nc**. For this class, the user may also turn off the slope limiting by modifying **limit_slopes**, but this is not recommended. 

Solver operator splitting
_________________________

Not all of our CDR solvers are compatible with our :ref:`Chap:time_stepper` implementations. The reason for this is that, internally in our code, we perform a segregation between CDR solvers. There is one category for the type of solvers that can perform a splitting of the convection and diffusion parts of :math:`\nabla\cdot(D\nabla n - \mathbf{v}n)`, and there is one category for solvers which cannot perform this decoupling. The Scharfetter-Gummel scheme, which computes the combined face flux as

.. math::

   F = v\frac{n_-\exp\left(v\Delta x/D\right) - n_+\exp\left(-v\Delta x/D\right)}{\exp\left(v\Delta x/D\right) - \exp\left(-v\Delta x/D\right)},

cannot be decoupled. This also means that the Scharfetter-Gummel scheme is *not* compatible with implementations of :ref:`Chap:time_stepper` that use operator splitting. On the other hand, our Godunov implementation uses slope-limited states for the advected state at the face, and *is* compatible with operator splitting methods. These solvers also have access to elliptic operators that compute :math:`\nabla\cdot\left(D\nabla n\right)`, and also the infrastructure that is used for implicit diffusion. 

In PlasmaC, new implementations of :doxy:`cdr_solver <classcdr__solver>` must implement the following routines:

.. code-block:: c++
		
  /*!
    @brief Compute div(nv - D*grad(n))
  */
  virtual void compute_divJ(EBAMRCellData& a_divJ, const EBAMRCellData& a_state, const Real a_extrap_dt) = 0;

  /*!
    @brief Compute advective derivative
  */
  virtual void compute_divF(EBAMRCellData& a_divF, const EBAMRCellData& a_state, const Real a_extrap_dt, const bool a_redist) = 0;

  /*!
    @brief Compute diffusion term
  */
  virtual void compute_divD(EBAMRCellData& a_diffusive_term, const EBAMRCellData& a_state) = 0;

The three functions above are used by the existing implementations of the CDR solvers. Here, *compute_divJ* should compute the full operator, whereas *compute_divF* and *compute_divD* allow segregation of advection and diffusion computations. For example, TGA implementations implement *compute_divJ* as

.. code-block:: c++

		void cdr_tga::compute_divJ(EBAMRCellData& a_divJ, const EBAMRCellData& a_state, const Real a_extrap_dt){
		   CH_TIME("cdr_tga::compute_divJ(divF, state)");
		   if(m_verbosity > 5){
		      pout() << m_name + "::compute_divJ(divF, state)" << endl;
		   }
		   
		   const int comp  = 0;
		   const int ncomp = 1;
		   
		   data_ops::set_value(a_divJ, 0.0);
		   
		   EBAMRCellData advective_term;
		   EBAMRCellData diffusion_term;
		   m_amr->allocate(advective_term, m_phase, ncomp);
		   if(this->is_diffusive()){
		      m_amr->allocate(diffusion_term, m_phase, ncomp);
		   }
		   
		   // Compute advective term
		   if(this->is_mobile()){
		      this->compute_divF(advective_term, a_state, 0.0, true);
		      data_ops::incr(a_divJ, advective_term, 1.0);
		   }
		   
		   // Add in diffusion term
		   if(this->is_diffusive()){
		      this->compute_divD(diffusion_term, a_state); // This already does refluxing. 
		      data_ops::incr(a_divJ, diffusion_term, -1.0);
		   }
		}

The internal workings of the *compute_divF* and *compute_divD* functions are more complicated and cannot be given in detail here. However, we mention that the *compute_divF* function computes the slope-limited face-centered states and the cell fluxes. The final update of that function also contains the post-processing that is involved for embedded boundaries (refluxing and redistribution). The elliptic part, *compute_divD*, essentially only contains a single update of a variable-coefficient Poisson equation.



