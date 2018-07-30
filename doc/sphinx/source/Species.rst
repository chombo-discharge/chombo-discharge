.. _Chap:species:

species
-------

The :ref:`Chap:species` is a lightweight class used to provide information into convection-diffusion-reaction solvers. This class is mostly used within :ref:`Chap:plasma_kinetics` in order to provide information on how to instantiate CDR solvers. :ref:`Chap:species` is abstract so that the user must implement

.. code-block:: c++

   /*!
    @brief Initial data. 
    @param[in] a_pos Position
    @param[in] a_time Time
  */
  virtual Real initial_data(const RealVect a_pos, const Real a_time) const = 0;

This function specifies the initial data of the species that is advected. For example, the following implementation sets the initial CDR density value to one:

.. code-block:: c++
		
  Real initial_data(const RealVect a_pos, const Real a_time) const {
     return 1.0;
  }

In addition to this, the user *must* provide information on the charge of the species, and whether or not it is mobile or diffusive. In addition, he should set the name of the species so that it can be identified in output files. In PlasmaC, this is done by setting the following four values in the constructor

.. code-block:: c++
		
  /*!
    @brief Species name
  */
  std::string m_name;
  
  /*!
    @brief Charge
  */
  int m_charge;

  /*!
    @brief Diffusive species or not
  */
  bool m_diffusive;

  /*!
    @brief Mobile species or not
  */
  bool m_mobile;


Usually, these are set through the constructor. The **m_charge** unit is in units of the elementary charge. For example, the following is a full implementation of an electron species:


.. code-block:: c++

		class electron : public species {
		  electron() {
		     m_name   = "electrons";
		     m_charge = -1;
		     m_diffusive = true;
		     m_mobile = true;
		  }

		  ~electron(){}

		  Real initial_data(const RealVect a_pos, const Real a_time) const {
		     return 1.0;
		  }
		};



The members **m_mobile** and **m_diffusive** are used for optimization in PlasmaC: If the user specifies that a species is immobile, PlasmaC will skip the advection computation. Note that **m_diffusive** and **m_mobile** override the specifications in :ref:`Chap:plasma_kinetics`. If the user provides a non-zero velocity through :ref:`Chap:plasma_kinetics` function *compute_cdr_velocities*, and sets **m_mobile** to **false**, the species velocity will be zero. Of course, the user will often want to provide additional input information to his species, for example by specifying a seed for the initial conditions. Please see :ref:`Chap:MiniApplications` for how to provide input parameters. 
