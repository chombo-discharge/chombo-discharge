.. _Chap:ItoKMC:

Îto-KMC plasma model
********************

.. warning::

   This section is not very well documented. The following featurse are definitely missing:

   * Use of hybrid models (PPC-fluid specifications)
   * How the physics time step is restricted

Underlying model
================

Plasma transport
----------------

The Îto-KMC model uses an Îto solver for some of the species, i.e. the particles drift and diffuse according to

.. math::

   d\mathbf{X} = \mathbf{V}dt + \sqrt{2D dt},

where :math:`\mathbf{X}` is the particle position and :math:`\mathbf{V}` and :math:`D` is the particle drift velocity and diffusion coefficients, respectively.
These are obtained by interpolating the macroscopic drift velocity and diffusion coefficients to the particle position.
Further details regarding the Îto method are given in :ref:`Chap:ItoSolver`.

Not all species in the Îto-KMC model need to be defined by particle solvers, as some species can be tracked by more conventional advection-diffusion-reaction solvers, see :ref:`Chap:CdrSolver` for discretization details.


Photoionization
---------------

Photoionization in the Îto-KMC model is done using discrete photons.
These are generated and advanced in every time step, and interact with the plasma through user-defined photo-reactions.
The underlying solver is the discrete Monte Carlo photon transport solver, see :ref:`Chap:MonteCarloRTE`.

Field interaction
-----------------

The Îto-KMC model uses an electrostatic approximation where the field is obtained by solving the Poisson equation for the potential.
See :ref:`Chap:Electrostatics` for further details.

Chemistry
---------

Kinetic Monte Carlo (KMC) is used within grid cells for resolving the plasma chemistry.
The algorithmic concepts are given in :ref:`Chap:KineticMonteCarlo`.

Algorithms
==========

Time stepping
-------------

The Îto-KMC model has adaptive time step controls and no fundamental CFL limitation.
Currently, the time step can be restricted through several means:

#. A physics-based time step, which is proportional to the non-critical time step from the KMC integrator. 
#. A particle advective, diffusive, or advective-diffusive CFL limitation, both upper and lower bounds.
#. A CFL-condition on the convection-diffusion-reaction solvers.
#. A limit that restricts the time step to a factor proportional to the dielectric relaxation time.
#. Hard limits, placing upper and lower bounds on the time step.

In addition, the user can specify the maximum permitted growth or reduction in the time step.

These limits are given by the following input variables:

.. literalinclude:: ../../../../Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepper.options
   :language: text
   :lines: 22-35

Particle placement
------------------

The KMC algorithm resolves the number of particles that are generated or lost within each grid cell, leaving substantial freedom in how one distributes new particles (or remove older ones). We currently support three methods for placing new particles:

#. Place all particles on the cell centroid.
#. Randomly distribute the new particles within the grid cell.
#. Compute an upstream position within the grid cell and randomly distribute the particles in the downstream region.

The methods each have their advantages and disadvantages.
Placing all new particles on the cell centroid has the advantage that there will be no spurious space charge effects arising in the grid cell.
Randomly distributing the new particles has the advantage that it generates secondary particles more evenly over the reaction region.
Both these methods (``centroid`` and ``random``) are, however, sources of numerical diffusion since the secondary particles are potentially placed in the wake of the primary particles (which is non-physical).
The downstream method circumvents this source of numerical diffusion by only placing secondary particles in the downstream region of some user-defined species (typically the electrons).
See :ref:`Chap:ItoKMCJSON` for instructions on how to assign the particle placement method.

Parallel diffusion
------------------

``ItoKMCGodunovStepper`` can limit diffusion against the electric field (or strictly speaking, the particle drift direction) by setting the flag ``ItoKMCGodunovStepper.limit_parallel_diffusion`` to ``true``.
The logic behind this functionality is that in a drift-diffusion approximation, and regardless of whether or one uses particles or fluids, electrons that back-diffuse against the electric field will rapidly lose their energy and can then no longer ionize the gas.
When limiting parallel diffusion, the diffusion hop of the particles is restricted to being along or transverse to the drift direction.
In practice, this leads to greatly enhanced stability in, especially in cathode sheaths.


Spatial filtering
-----------------

It is possible to apply filtering of the space-charge density prior to the field advance in order to reduce the impact of discrete particle noise.
Filters are applied as follows:

.. math::

   \rho_i = \alpha\rho_i + \frac{1-\alpha}{2}\left(\rho_{i+s} + \rho_{i-s}\right).

where :math:`\alpha` is a filtering factor and :math:`s` a stride.
Users can apply this filtering by adjusting the following input options:

.. code-block:: text

   ItoKMCGodunovStepper.rho_filter_num        = 0    # Number of filterings for the space-density
   ItoKMCGodunovStepper.rho_filter_max_stride = 1    # Maximum stride for filter
   ItoKMCGodunovStepper.rho_filter_alpha      = 0.5  # Filtering factor (0.5 is a bilinear filter)		

.. warning::

   Spatial filtering of the space-charge density is a work-in-progress which may yield unexpected results.
   Bugs may or may not be present. 
   Users should exercise caution when using this feature.



Reaction network
----------------

Particle management
-------------------

0D chemistry
============

The user input interface to the Îto-KMC model consists of a zero-dimensional plasma kinetics interface called ``ItoKMCPhysics``.
This interface consists of the following main functionalities:

#. User-defined species, and which solver types (Îto or CDR) to use when tracking them.
#. User-defined initial particles, reaction kinetics, and photoionization processes.

.. _Chap:ItoKMCPhysics:

ItoKMCPhysics
-------------

The complete C++ interface specification is given below.
Because the interface is fairly extensive, ``chombo-discharge`` also supplies a JSON-based implementation called ``ItoKMCJSON`` (see :ref:`Chap:ItoKMCJSON`) for defining these things through file input.
The difference between ``ItoKMCPhysics`` and its implementation ``ItoKMCJSON`` is that the JSON-based implementation class only implements a subset of potential features supported by ``ItoKMCPhysics``.

.. raw:: html

   <details>
   <summary><a>Show/hide full 0D interface</a></summary>

.. literalinclude:: ../../../../Physics/ItoKMC/CD_ItoKMCPhysics.H
   :language: c++   

.. raw:: html

   </details><br>   

Species definitions
-------------------

Species are defined either as input species for CDR solvers (see :ref:`Chap:CdrSolver`) or Îto solvers (see :ref:`Chap:ItoSolver`).
It is sufficient to populate the ``ItoKMCPhysics`` species vectors ``m_cdrSpecies`` and ``m_itoSpecies`` if one only wants to initialize the solvers.
See :ref:`Chap:ItoKMCPhysics` for their C++ definition.
In addition to actually populating the vectors, users will typically also include initial conditions for the species.
The interfaces permit initial particles and density functions for both both solver types.
Additionally, one must define all species associated with radiative transfer solvers.

.. _Chap:ItoKMCPlasmaReaction:

Plasma reactions
----------------

Plasma reactions in the Îto-KMC model are represented stoichiometrically as

.. math::

   A + B + \ldots \xrightarrow{c} C + D + \ldots

where :math:`c` is the KMC rate coefficient, i.e. *not the rate coefficient from the reaction rate equation*.
An arbitrary number of reactions is supported, but currently the reaction mechanisms are limited to:

* The left-hand side must consist only of species that are tracked by a CDR or Îto solver.
* The right-hand side can consist of species tracked by a CDR solver, an Îto solver, or a radiative transfer solver.

These rules apply only to the internal C++ interface; implementations of this interface will typically also allow species defined as ''background species'', i.e. species that are not tracked by a solver.
In that case, the interface implementation must absorb the background species contribution into the rate constant.

In the KMC algorithm, one operates with chemical propensities rather than rate coefficients.
For example, the chemical propensity :math:`a_1` for a reaction :math:`A + B \xrightarrow{c_1}\varnothing` is

.. math::

   a_1 = c_1X_AX_B,

and from the macroscopic limit :math:`X_A \gg 1`, :math:`X_B \gg 1` one may in fact also derive that :math:`c_1 = k_1/\Delta V` where :math:`k_1` is the rate coefficient from the reaction rate equation, i.e. where :math:`\partial_t n_A = -k_1n_An_B`.
Note that if the reaction was :math:`A + A\xrightarrow{c_2}\varnothing` then the propensity is

.. math::

   a_2 = \frac{1}{2}c_2X_A\left(X_A-1\right)

since there are :math:`\frac{1}{2}X_A\left(X_A-1\right)` distinct pairs of particles of type :math:`A`.
Likewise, the fluid rate coefficient would be :math:`k_2 = c_2\Delta V/2`.

The distinction between KMC and fluid rates is an important one; the reaction representation used in the Îto-KMC model only operates with the KMC rates :math:`c_j`, and it is us to the user to ensure that these are consistent with the fluid limit.
Internally, these reactions are implemented through the dual state KMC implementation, see :ref:`Chap:KineticMonteCarlo`.
During the reaction advance the user only needs to update the :math:`c_j` coefficients (typically done via an interface implementation); the calculation of the propensity is automatic and follows the standard KMC rules (e.g., the KMC solver accounts for the number of distinct pairs of particles).
This must be done in the routine ``updateReactionRates(...)``, see :ref:`Chap:ItoKMCPhysics` for the complete specification.

Photoionization
---------------

Photo-reactions are also represented stoichiometrically as

.. math::

   \gamma \xrightarrow{\xi} A + B + \ldots,

where :math:`\xi` is the photo-reaction efficiency, i.e. the probability that the photon :math:`\gamma` causes a reaction when it is absorbed on the mesh.
Currently, photons can only be generated through *plasma reactions*, see :ref:`Chap:ItoKMCPlasmaReaction`, and we do not permit reactions between the photon and plasma species.

.. important::

   We currently only support constant photoionization probabilities :math:`\xi`.

``ItoKMCPhysics`` adds some flexibility when dealing with photo-reactions as it permits pre-evaluation of photoionization probabilities.
That is, when generating photons through reactions of the type :math:`A + B \rightarrow \gamma`, one may scale the reaction by the photoionization probability :math:`\xi` so that only ionizing photons are generated and then write the photo-reaction as :math:`\gamma \xrightarrow{\xi=1} A + B + \ldots`.
However, this process becomes more complicated when dealing with multiple photo-reaction pathways, e.g.,

.. math::

   \gamma &\xrightarrow{\xi_1} A + B, \\
   \gamma &\xrightarrow{\xi_2} A + B, \\

where :math:`\xi_1` and :math:`\xi_2` are the probabilities that the photon :math:`\gamma` triggers the reaction.
If only a single physical photon is absorbed, then one must stochastically determine which pathway is triggered.
There are two ways of doing this:

#. Pre-scale the photon generation by :math:`\xi_1 + \xi_2`, and then determine the pathway through the relative probabilities :math:`\xi_1/\left(\xi_1+\xi_2\right)` and :math:`\xi_2/\left(\xi_1+\xi_2\right)`.
#. Do not scale the photon generation reaction and determine the pathway through the absolute probabilities :math:`\xi_1` and :math:`\xi_2`.
   This can imply a large computational cost since one will have to track all photons that are generated.

The former method is normally the preferred way as it can lead to reductions in computational complexity and particle noise, but for flexibility ``ItoKMCPhysics`` supports both of these.
The latter method can be of relevance if users wants precise descriptions of photons that trigger both photoionizing reactions and surface reactions (e.g., secondary electron emission).
   
An alternative is to split the photon type :math:`\gamma` into volumetrically absorbed and surface absorbed photon species :math:`\gamma_v` and :math:`\gamma_s`.
In this case one may pre-scale the photon-generating reactions by the photo-reaction probabilities :math:`\xi` (for volume absorption) and :math:`\zeta` (for surface absorption) as follows:

.. math::

   A + B &\xrightarrow{k\rightarrow k\xi}\gamma_v, \\ 
   A + B &\xrightarrow{k\rightarrow k\zeta}\gamma_s.

Volumetric and surface absorption is then treated independently

.. math::

   \gamma_v &\xrightarrow{\left(\xi=1\right)}, \ldots\\
   \gamma_s &\xrightarrow{\left(\zeta=1\right)} \ldots.

This type of pre-evaluation of the photo-reaction pathways is sensible in a statistical sense, but loses meaning if only a single photon is involved.

Surface reactions
-----------------

.. warning::

   Surface reactions are supported by :ref:`Chap:ItoKMCPhysics` but not implemented in the JSON interface (yet).

Transport coefficients
----------------------

Species mobilities and diffusion coefficients should be computed just as they are done in the fluid approximation.
The ``ItoKMCPhysics`` interface requires implementations of two functions that define the coefficients as functions of :math:`\mu = \mu\left(t,\mathbf{x}, \mathbf{E}\right)`,
and likewise for the diffusion coefficients.
Note that these functions should return the *fluid coefficients*.

.. important::

   There is currently no support for computing :math:`\mu` as a function of the species densities (e.g., the electron density), but this only requires modest extensions of the Îto-KMC module.

Time step calculation
=====================

.. warning::

   This section must be written, and also include the downsides of using :math:`\Delta t = X/|\sum\mu|` in its pure form. 

.. _Chap:ItoKMCJSON:

JSON interface
==============

The JSON-based chemistry interface (called ``ItoKMCJSON``) simplifies the definition of the plasma kinetics and initial conditions by permitting specifications through a JSON file.
In addition, the JSON specification will permit the definition of background species that are not tracked as solvers (but that simply exist as a density function).
The JSON interface thus provides a simple entry point for users that have no interest in defining the chemistry from scratch.

Several mandatory fields are required when specifying the plasma kinetics, such as:

* Gas pressure, temperature, and number density.
* Background species (if any).
* Townsend coefficients.
* Plasma species and their initial conditions, which solver to use when tracking them.
* Photon species and their absorption lengths.
* Plasma reactions.
* Photoionization reactions. 

These specifications follow a pre-defined JSON format which is discussed in detail below.

.. tip::

   While JSON files do not formally support comments, the JSON submodule used by ``chombo-discharge`` allows them.

Gas law
-------

The JSON gas law represents a loose definition of the gas pressure, temperature, and *number density*.
Multiple gas laws can be defined in the JSON file, and the user can then select which one to use.
This is defined within a JSON entry ``gas`` and a subentry ``law``.
In the ``gas/law`` field one must specify an ``id`` field which indicates which gas law to use.
One may also set an optional field ``plot`` to true/false if one wants to include the pressure, temperature, and density in the HDF5 output files.

.. important::

   Gas parameters are specified in SI units.

Each user-defined gas law exists as a separate entry in the ``gas/law`` field, and the minimum requirements for these are that they contain a ``type`` specifier.
Currently, the supported types are

* ``ideal`` (for constant ideal gas)
  
An example JSON specification is given below.

.. raw:: html

   <details>
   <summary><a>Example gas law specification</a></summary>

.. code-block:: json

   {
      "gas" : {
         // Specify which gas law we use. The user can define multiple gas laws and then later specify which one to use.
	 "law" : {
	    "id" : "my_ideal_gas", // Specify which gas law we use.
	    "plot" : true,         // Turn on/off plotting.
	    "my_ideal_gas" : {
		// Definition for an ideal gas law. Temperature is always in Kelvin the pressure is in Pascal. The neutral density
		// is derived using an ideal gas law.
		"type" : "ideal",
		"temperature" : 300,
		"pressure" : 1E5
	    }
	 }
      }
   }

.. raw:: html

   </details><br>


Ideal gas
_________

When specifying that the gas law ``type`` is ``ideal``, the user must further specify the temperature and pressure of the gas.

.. raw:: html

   <details>
   <summary><a>Example ideal gas specification</a></summary>

.. code-block:: json

   {
      "gas" : {
	 "law" : {
	    "id" : "my_ideal_gas", // Specify which gas law we use.
	    "my_ideal_gas" : {
		"type" : "ideal",    // Ideal gas law type
		"temperature" : 300, // Temperature must be in Kelvin
		"pressure" : 1.0     // Pressure must be in bar
	    },
	 }
      }
   }


.. raw:: html

   </details><br>


..
   Table vs height
   _______________

   .. warning::

      Not yet implemented.

Background species
------------------

Background species :math:`i` are defined as continuous background densities :math:`N_i` given by

.. math::

   N_i\left(\mathbf{x}\right) = \chi_i\left(\mathbf{x}\right)N\left(\mathbf{x}\right),

where :math:`\chi_i\left(\mathbf{x}\right)` is the molar fraction of species :math:`i` (typically, but not necessarily, a constant).

When specifying background species, the user must include an array of background species in the ``gas`` JSON field.
For example,

.. code-block:: json

   {
      "gas" : {
	 "background species" : [
	    {
	       // First species goes here
	    },
	    {
	       // Second species goes here
	    }	    
	 ]
      }
   }

The order of appearance of the background species does not matter. 

Name
____

The species name is specified by including an ``id`` field, e.g.

.. code-block:: json
   :emphasize-lines: 5
      
   {
      "gas" : {
	 "background species" : [
	    {
	       "id": "O2"
	    }
	 ]
      }
   }

Molar fraction
______________

The molar fraction can be specified in various forms by including a field ``molar fraction``.
This field also requires the user to specify the *type* of molar fraction.
In principle, the molar fraction can have a positional dependency.

.. warning::

   There's no internal renormalization of the molar fractions input by the user.
   Internal inconsitencies will occur if the user supplies inputs molar fractions that sum to a number different than one.

   Various checks will be enabled if the user compiles in debug-mode (see :ref:`Chap:Installation`), but these checks are not guaranteed to catch all cases.

Constant
^^^^^^^^

To specify a constant molar fraction, set the ``type`` specifier to constant and specify the value.
An example JSON specification is

.. code-block:: json
   :emphasize-lines: 6-9
      
   {
      "gas" : {
	 "background species" : [
	    {
	       "id": "O2"             // Species name
	       "molar fraction": {    // Molar fraction specification
   	          "type": "constant", // Constant molar fraction
	          "value": 0.2        // Molar fraction value
	       }
	    }
	 ]
      }
   }

Tabulated versus height
^^^^^^^^^^^^^^^^^^^^^^^

The molar fraction can be set as a tabulated value versus one of the Cartesian coordinate axis by setting the ``type`` specifier to ``table vs height``.
The input data should be tabulated in column form, e.g.

.. code-block:: text

   # height       molar fraction
   0              0.1
   1              0.1 
   2              0.1

The file parser (see :ref:`Chap:LookupTable`) will ignore the header file if it starts with a hashtag (#).
Various other inputs are then also required:

* ``file`` File name containing the height vs. molar fraction data (required).
* ``axis`` Cartesian axis to associate with the height coordinate (required).
* ``height column`` Column in data file containing the height (optional, defaults to 0).
* ``molar fraction column`` Column in data file containing the molar fraction (optional, defaults to 1).
* ``height scale`` Scaling of height column (optional).
* ``fraction scale`` Scaling (optional).
* ``min height`` Truncate data to minimum height (optional, applied after scaling).
* ``max height`` Truncate data to maximum height (optional, applied after scaling).
* ``num points`` Number of points to keep in table representation (optional, defaults to 1000).
* ``spacing`` Spacing table representation (optional, defaults to "linear".
* ``dump`` Option for dumping tabulated data to file.

An example JSON specification is

.. code-block:: json
   :emphasize-lines: 6-19
      
   {
      "gas" : {
	 "background species" : [
	    {
	       "id": "O2"                           // Species name
	       "molar fraction": {                  // Molar fraction specification
   	          "type": "table vs height",        // Constant molar fraction
		  "file" : "O2_molar_fraction.dat", // File name containing molar fraction
		  "dump" : "debug_O2_fraction.dat", // Optional debugging hook for dumping table to file
		  "axis" : "y",                     // Axis which represents the "height"
		  "height column" : 0,              // Optional specification of column containing the height data (defaults to 0)
		  "molar fraction column" : 1,      // Optional specification of column containing the height data (defaults to 1)
		  "height scale" : 1.0,             // Optional scaling of height column
		  "fraction scale" : 1.0,           // Optional scaling of molar fraction column
		  "min height" : 0.0,               // Optional truncation of minimum height kept in internal table (occurs after scaling)
		  "max height" : 2.0,               // Optional truncation of maximum height kept in internal table (occurs after scaling)
		  "num points" : 100,               // Optional specification of number of data points kept in internal table (defaults to 1000)
		  "spacing" : "linear"              // Optional specification of table representation. Can be 'linear' or 'exponential' (defaults to linear)
	       }
	    }
	 ]
      }
   }

Plotting
________

Background species densities can be included in HDF5 plots by including an optional field ``plot``.
For example

.. code-block:: json
   :emphasize-lines: 10
      
   {
      "gas" : {
	 "background species" : [
	    {
	       "id": "O2"             // Species name
	       "molar fraction": {    // Molar fraction specification
   	          "type": "constant", // Constant molar fraction
	          "value": 0.2        // Molar fraction value
	       },
	       "plot": true           // Plot number density in HDF5 or not
	    }
	 ]
      }
   }

Townsend coefficients
---------------------

Townsend ionization and attachment coefficients :math:`\alpha` and :math:`\eta` must be specified, and have two usages:

#. Flagging cells for refinement (users can override this).
#. Potential usage in plasma reactions (which requires particular care, see :ref:`Chap:ItoKMCWarnings`. 

These are specified by including JSON entries ``alpha`` and ``eta``, e.g.

.. code-block:: json
	
   {
       "alpha": {
          // alpha-coefficient specification goes here
       },
       "eta" : {
          // eta-coefficient specification goes here
       }
   }

There are various way of specifying these, as discussed below:

Automatic
_________

In auto-mode the Townsend coefficients for the electrons is automatically derived from the user-specified list of reactions.
E.g., the Townsend ionization coefficient is computed as

.. math::

   \alpha = \frac{\sum k}{\mu E},

where :math:`\sum k` is the sum of all ionizing reactions, also incorporating the neutral density.
The advantage of this approach is that one may choose to discard some of the reactions from e.g. BOLSIG+ output but without recomputing the Townsend coefficients.
   
To automatically set Townsend coefficients, set ``type`` to auto and specify the species that is involved, e.g.

.. code-block:: json
		
   {
      "alpha": {
         "type": "auto",
	 "species": "e"
      }
   }

The advantage of using auto-mode for the coefficients is that one automatically ensures consistency between the user-specified list of reactions and the 
   

Constant
________

To set a constant Townsend coefficient, set ``type`` to constant and then specify the value, e.g.

.. code-block:: json
   :emphasize-lines: 3-4
      
   {
       "alpha": {
          "type": "constant",
	  "value": 1E5
       }
   }


Tabulated vs :math:`E/N`
________________________
    
To set the coefficient as functions :math:`f = f\left(E/N\right)`, set the ``type`` specifier to ``table vs E/N`` and specify the following fields:

* ``file`` For specifying the file containing the input data, which must be organized as column data, e.g.

  .. code-block:: text

     # E/N   alpha/N
     0       1E5
     100     1E5
     200     1E5

* ``header`` For specifying where in the file one starts reading data.
  This is an optional argument intended for use with BOLSIG+ output data where the user specifies that the data is contained in the lines below the specified header.
* ``E/N column`` For setting which column in the data file contains the values of :math:`E/N` (optional, defaults to 0).
* ``alpha/N column`` For setting which column in the data file contains the values of :math:`\alpha/N`.
  Note that this should be replaced by ``eta/N column`` when parsing the attachment coefficient.
  This is an optional argument that defaults to 1.
* ``scale E/N`` Optional scaling of the column containing the :math:`E/N` data.
* ``scale alpha/N`` Optional scaling of the column containing the :math:`\alpha/N` data.
* ``min E/N`` Optional truncation of the table representation of the data (applied after scaling).
* ``max E/N`` Optional truncation of the table representation of the data (applied after scaling).
* ``num points`` Number of data points in internal table representation (optional, defaults to 1000).
* ``spacing`` Spacing in internal table representation. Can be ``linear`` or ``exponential``.
  This is an optional argument that defaults to ``exponential``.
* ``dump`` A string specifier for writing the table representation of :math:`E/N, \alpha/N` to file.

An example JSON specification that uses a BOLSIG+ output file for parsing the data is

.. code-block:: json
   :emphasize-lines: 3-14
      
   {
       "alpha": {
          "type" : "table vs E/N",                                   // Specify that we read in alpha/N vs E/N
  	  "file" : "bolsig_air.dat",                                 // File containing the data
	  "dump" : "debug_alpha.dat",                                // Optional dump of internalized table to file. Useful for debugging.	
	  "header" : "E/N (Td)\tTownsend ioniz. coef. alpha/N (m2)", // Optional argument. Contains line immediately preceding the data to be read.
	  "E/N column" : 0,                                          // Optional specification of column containing E/N (defaults to 0)
	  "alpha/N column" : 1,                                      // Optional specification of column containing alpha/N (defaults to 1)
	  "min E/N" : 1.0,                                           // Optional truncation of minium E/N kept when resampling the table (occurs after scaling)
	  "max E/N" : 1000.0,                                        // Optional truncation of maximum E/N kept when resampling the table (happens after scaling)
	  "num points" : 1000,                                       // Optional number of points kept when resamplnig the table (defaults to 1000)
	  "spacing" : "exponential",                                 // Optional spcification of table representation. Defaults to 'exponential' but can also be 'linear'
	  "scale E/N" : 1.0,                                         // Optional scaling of the column containing E/N
	  "scale alpha/N" : 1.0                                      // Optional scaling       
       }
   }

Plotting
________

To include the Townsend coefficients as mesh variables in HDF5 files, include the ``plot`` specifier, e.g.

.. code-block:: json
   :emphasize-lines: 5
      
   {
      "alpha": {
         "type": "auto",
	 "species": "e",
	 "plot": true
      }
   }

Plasma species
--------------

Plasma species are species that are tracked using either a CDR or Îto solver.
The user must specify the following information:

* An ID/name for the species.
* The species' charge number.
* The solver type, i.e. whether or not it is tracked by a particle or fluid solver. 
* Whether or not the species is mobile/diffusive.
  If the species is mobile/diffusive, the mobility/diffusion coefficient must also be specified.
* Optionally, the species temperature.
  If not specified, the temperature will be set to the background gas temperature. 
* Initial particles for the solvers. 


Basic definition
________________

Species are defined by an array of entries in a ``plasma species`` JSON field.
Each species *must* specify the following fields

* ``id`` (string) For setting the species name.
* ``Z`` (integer) For setting the species' charge number.
* ``solver`` (string) For specifying whether or not the species is tracked by a CDR or Îto solver.
  Acceptable values are *cdr* and *ito*.
* ``mobile`` (boolean) For specifying whether or not the species is mobile.
* ``mobile`` (boolean) For specifying whether or not the species is diffusive.

An example specification is

.. code-block:: json

    "plasma species" :
    [
	// List of plasma species that are tracked. This is an array of species
	// with various identifiers, some of which are always required (id, Z, type, mobile, diffusive) and
	// others which are secondary requirements.
	{
            "id": "e",          // Species ID
	    "Z" : -1,           // Charge number
	    "solver" : "ito",   // Solver type. Either 'ito' or 'cdr'
	    "mobile" : true,    // Mobile or not
	    "diffusive" : true  // Diffusive or not
	},
	{
	    // Definition of O2+ plasma species.
            "id": "O2+",            // Species ID
	    "Z" : 1,                // Charge number
	    "solver" : "cdr",       // CDR solver. 
	    "mobile" : false,       // Not mobile.
	    "diffusive" : false     // Not diffusive
	}
    ]

Note that the order of appearance of the various species is irrelevant.
However, if a species is specified as mobile/diffusive, the user *must* also specify the corresponding transport coefficients. 

Mobility and diffusion coefficients
___________________________________

Mobility and diffusion coefficients are specified by including fields ``mobility`` and ``diffusion`` that further specifies how the transport coefficients are calculated.
Note that the input to these fields should be *fluid transport coefficients*.
These fields can be specified in various forms, as shown below.

Constant
^^^^^^^^

To set a constant coefficient, set the ``type`` specifier to constant and then assign the value.
For example,

.. code-block:: json
   :emphasize-lines: 10-11

                     
    "plasma species" :
    [
	{
            "id": "e",             // Species ID
	    "Z" : -1,              // Charge number
	    "solver" : "ito",      // Solver type. Either 'ito' or 'cdr'
	    "mobile" : true,       // Mobile or not
	    "diffusive" : true     // Diffusive or not,
	    "mobility": {
	       "type": "constant", // Set constant mobility
	       "value": 0.02
	    },
	    "diffusion": {
	       "type": "constant", // Set constant diffusion coefficient
	       "value": 2E-4
	    }
	}
    ]

Constant :math:`*N`
^^^^^^^^^^^^^^^^^^^^^^

To set a coefficient that is constant vs :math:`N`, set the ``type`` specifier to ``constant mu*N`` or ``constant D*N`` and then assign the value.
For example,

.. code-block:: json
   :emphasize-lines: 10-11

    "plasma species" :
    [
	{
            "id": "e",			// Species ID
	    "Z" : -1,			// Charge number
	    "solver" : "ito",		// Solver type. Either 'ito' or 'cdr'
	    "mobile" : true,		// Mobile or not
	    "diffusive" : true		// Diffusive or not,
	    "mobility": {
	       "type": "constant mu*N", // Set mu*N to a constant
	       "value": 1E24
	    },
	    "diffusion": {
	       "type": "constant D*N",  // Set D*N to a constant
	       "value": 5E24
	    }
	}
    ]    

Table vs :math:`E/N`
^^^^^^^^^^^^^^^^^^^^

To set the transport coefficients as functions :math:`f = f\left(E/N\right)`, set the ``type`` specifier to ``table vs E/N`` and specify the following fields:

* ``file`` For specifying the file containing the input data, which must be organized as column data, e.g.

  .. code-block:: text

     # E/N   mu/N
     0       1E5
     100     1E5
     200     1E5

* ``header`` For specifying where in the file one starts reading data.
  This is an optional argument intended for use with BOLSIG+ output data where the user specifies that the data is contained in the lines below the specified header.
* ``E/N column`` For setting which column in the data file contains the values of :math:`E/N` (optional, defaults to 0).
* ``mu*N column`` For setting which column in the data file contains the values of :math:`\mu N` (or alternatively :math:`DN` for the diffusion coefficient).
  This is an optional argument that defaults to 1.
* ``E/N scale`` Optional scaling of the column containing the :math:`E/N` data.
* ``mu*N scale`` Optional scaling of the column containing the :math:`\mu N` data (or alternatively :math:`DN` for the diffusion coefficient).
* ``min E/N`` Optional truncation of the table representation of the data (applied after scaling).
* ``max E/N`` Optional truncation of the table representation of the data (applied after scaling).
* ``num points`` Number of data points in internal table representation (optional, defaults to 1000).
* ``spacing`` Spacing in internal table representation. Can be ``linear`` or ``exponential``.
  This is an optional argument that defaults to ``exponential``.
* ``dump`` A string specifier for writing the table representation of :math:`E/N, \mu N` to file.

An example JSON specification that uses a BOLSIG+ output file for parsing the data for the mobility is

.. code-block:: json
   :emphasize-lines: 9-22

    "plasma species" :
    [
       {
          "id": "e",                                       // Species ID
          "Z" : -1,                                        // Charge number
          "solver" : "ito",                                // Solver type. 
          "mobile" : true,                                 // Mobile
          "diffusive" : false,                             // Not diffusive
          "mobility" : {
	     "type" : "table vs E/N",                      // Specification of tabulated mobility lookup method
	     "file" : "bolsig_air.dat",                    // File containg the mobility data
	     "dump" : "debug_mobility.dat",                // Optional argument for dumping table to disk (useful for debugging)		
	     "header" : "E/N (Td)\tMobility *N (1/m/V/s)", // Line immediately preceding the colum data
	     "E/N column" : 0,                             // Column containing E/N
	     "mu*N column" : 1,                            // Column containing mu*N
	     "min E/N" : 1,                                // Minimum E/N kept when resampling table
	     "max E/N" : 2E6,                              // Maximum E/N kept when resampling table
	     "points" : 1000,                              // Number of grid points kept when resampling the table
	     "spacing" : "exponential",                    // Grid point spacing when resampling the table
	     "E/N scale" : 1.0,                            // Optional argument for scaling mobility coefficient
	     "mu*N scale" : 1.0                            // Optional argument for scaling the E/N column
	  }
       }
    ]

.. tip::

   The parser for the diffusion coefficient is analogous; simply replace ``mu*N`` by ``D*N``.
    

Temperature
___________

The species temperature is only parametrically attached to the species (i.e., not solved for), and can be specified by the user.
Note that the temperature is mostly relevant for transport coefficients (e.g., reaction rates).
The temperature can be specified through the ``temperature`` field for each species, and various forms of specifying this is available.

.. important::
   
   If the specifies temperature is *not* specified by the user, it will be set equal to the background gas temperature.

Background gas
^^^^^^^^^^^^^^

To set the temperature equal to the background gas temperature, set the ``type`` specifier field to ``gas``.
For example:

.. code-block:: json
   :emphasize-lines: 7-9

    "plasma species" :
    [
       {
          "id": "O2+",      // Species ID
          "Z" : 1,          // Charge number
          "solver" : "ito", // Solver type. 
	  "temperature": {  // Specify temperature
	     "type": "gas"  // Temperature same as gas
	  }
       }
    ]

Constant
^^^^^^^^

To set a constant temperature, set the ``type`` specifier to ``constant`` and then specify the temperature.
For example:

.. code-block:: json
   :emphasize-lines: 7-10

    "plasma species" :
    [
       {
          "id": "O2+",           // Species ID
          "Z" : 1,               // Charge number
          "solver" : "ito",      // Solver type. 
	  "temperature": {       // Specify temperature
	     "type": "constant", // Constant temperature
	     "value": 300        // Temperature in Kelvin
	  }
       }
    ]

Table vs :math:`E/N`
^^^^^^^^^^^^^^^^^^^^

To set the species temperature as a function :math:`T = T\left(E/N\right)`, set the ``type`` specifier to ``table vs E/N`` and specify the following fields:

* ``file`` For specifying the file containing the input data, which must be organized as column data *versus the mean energy*, e.g.

  .. code-block:: text

     # E/N   energy (eV)
     0       1
     100     2
     200     3

* ``header`` For specifying where in the file one starts reading data.
  This is an optional argument intended for use with BOLSIG+ output data where the user specifies that the data is contained in the lines below the specified header.
* ``E/N column`` For setting which column in the data file contains the values of :math:`E/N` (optional, defaults to 0).
* ``eV column`` For setting which column in the data file contains the energy vs :math:`E/N`.
  This is an optional argument that defaults to 1.
* ``E/N scale`` Optional scaling of the column containing the :math:`E/N` data.
* ``eV scale`` Optional scaling of the column containing the energy (in eV).
* ``min E/N`` Optional truncation of the table representation of the data (applied after scaling).
* ``max E/N`` Optional truncation of the table representation of the data (applied after scaling).
* ``num points`` Number of data points in internal table representation (optional, defaults to 1000).
* ``spacing`` Spacing in internal table representation. Can be ``linear`` or ``exponential``.
  This is an optional argument that defaults to ``exponential``.
* ``dump`` A string specifier for writing the table representation of :math:`E/N, eV` to file.

An example JSON specification that uses a BOLSIG+ output file for parsing the data for the mobility is

.. code-block:: json
   :emphasize-lines: 9-22

    "plasma species" :
    [
       {
          "id": "O2+",                                // Species ID
          "Z" : 1,                                    // Charge number
          "solver" : "ito",                           // Solver type. 
          "mobile" : true,                            // Not mobile
          "diffusive" : false,                        // Not diffusive
	  "temperature": {                            // Specification of temperature
	     "type": "table vs E/N",                  // Tabulated
	     "file": "bolsig_air.dat",                // File name
	     "dump": "debug_temperature.dat",         // Dump to file
	     "header" : "E/N (Td)\tMean energy (eV)", // Header preceding data
	     "E/N column" : 0,                        // Column containing E/N
	     "eV column" : 1,                         // Column containing the energy
	     "min E/N" : 10,                          // Truncation of table
	     "max E/N" : 2E6,                         // Truncation of table
	     "E/N scale": 1.0,                        // Scaling of input data
	     "eV scale": 1.0,                         // Scaling of input data
	     "spacing" : "exponential",               // Table spacing
	     "points" : 1000                          // Number of points in table
	  }
       }
    ]

Initial particles
_________________

Initial particles for species are added by including an ``initial particles`` field.
The field consists of an array of JSON entries, where each array entry represents various ways of adding particles.

.. important::

   The ``initial particles`` field is *incrementing*, each new entry will add additional particles.

For example, to add two particles with particle weights of 1, one may specify

.. code-block:: json
   :emphasize-lines: 9-22
		     
    "plasma species" :
    [
       {
          "id": "e",                 // Species ID
          "Z" : 1,                   // Charge number
          "solver" : "ito",          // Solver type. 
          "mobile" : false,          // Not mobile
          "diffusive" : false,       // Not diffusive
	  "initial particles": [     // Initial particles
	     {
	        "single particle": { // Single physical particle at (0,0,0)
		   "position": [0, 0, 0],
		   "weight": 1.0
		}
	     },
	     {
	        "single particle": { // Single physical particle at (0,1,0)
		   "position": [0, 1, 0],
		   "weight": 1.0
		}
	     }
	  ]
       }
    ]

The various supported ways of additional initial particles are discussed below.

.. important::

   If a solver is specified as a CDR solver, the initial particles are deposited as densities on the mesh using an NGP scheme.

Single particle
^^^^^^^^^^^^^^^

Single particles are placed by including a ``single particle`` JSON entry.
The user must specify

* ``position`` The particle position.
* ``weight`` The particle weight.
  
For example:

.. code-block:: json
   :emphasize-lines: 11-14
		     
    "plasma species" :
    [
       {
          "id": "e",                 // Species ID
          "Z" : 1,                   // Charge number
          "solver" : "ito",          // Solver type. 
          "mobile" : false,          // Not mobile
          "diffusive" : false,       // Not diffusive
	  "initial particles": [     // Initial particles
	     {
	        "single particle": { // Single physical particle at (0,0,0)
		   "position": [0, 0, 0],
		   "weight": 1.0
		}
	     }
	  ]
       }
    ]

Uniform distribution
^^^^^^^^^^^^^^^^^^^^

Particles can be randomly drawn from a uniform distribution (a rectangular box in 2D/3D) by including a ``uniform distribution``.
The user must specify

* ``low corner`` The lower left corner of the box.
* ``high corner`` The lower left corner of the box.
* ``num particles`` Number of computational particles that are drawn.
* ``weight`` Particle weights.

For example:

.. code-block:: json
   :emphasize-lines: 11-16
		     
    "plasma species" :
    [
       {
          "id": "e",                                  // Species ID
          "Z" : 1,                                    // Charge number
          "solver" : "ito",                           // Solver type. 
          "mobile" : false,                           // Not mobile
          "diffusive" : false,                        // Not diffusive
	  "initial particles": [                      // Initial particles
	     {
	        "uniform distribution" : {            // Particles inside a box
		   "low corner" : [ -0.04, 0, 0 ],    // Lower-left physical corner of sampling region
		   "high corner" : [ 0.04, 0.04, 0 ], // Upper-right physical corner of sampling region
		   "num particles": 1000,             // Number of computational particles
		   "weight" : 1                       // Particle weights
		} 	     
	     }
	  ]
       }
    ]

Sphere distribution
^^^^^^^^^^^^^^^^^^^

Particles can be randomly drawn inside a circle (sphere in 3D)  by including a ``sphere distribution``.
The user must specify

* ``center`` The sphere center.
* ``radius`` The sphere radius.
* ``num particles`` Number of computational particles that are drawn.
* ``weight`` Particle weights.

For example:

.. code-block:: json
   :emphasize-lines: 11-16
		     
    "plasma species" :
    [
       {
          "id": "e",                            // Species ID
          "Z" : 1,                              // Charge number
          "solver" : "ito",                     // Solver type. 
          "mobile" : false,                     // Not mobile
          "diffusive" : false,                  // Not diffusive
	  "initial particles": [                // Initial particles
	     {
	        "sphere distribution" : {       // Particles inside sphere
		   "center" : [ 0, 7.5E-3, 0 ], // Sphere center
		   "radius" : 1E-3,             // Sphere radius		
		   "num particles": 1000,       // Number of computational particles
		   "weight" : 1                 // Particle weights
		} 	     
	     }
	  ]
       }
    ]

Gaussian distribution
^^^^^^^^^^^^^^^^^^^^^

Particles can be randomly drawn from a Gaussian distribution by including a ``gaussian distribution``.
The user must specify

* ``center`` The sphere center.
* ``radius`` The sphere radius.
* ``num particles`` Number of computational particles that are drawn.
* ``weight`` Particle weights.

For example:

.. code-block:: json
   :emphasize-lines: 11-16
		     
    "plasma species" :
    [
       {
          "id": "e",                            // Species ID
          "Z" : 1,                              // Charge number
          "solver" : "ito",                     // Solver type. 
          "mobile" : false,                     // Not mobile
          "diffusive" : false,                  // Not diffusive
	  "initial particles": [                // Initial particles
	     {
	        "gaussian distribution" : {     // Particles drawn from a gaussian 
		   "center" : [ 0, 7.5E-3, 0 ], // Center
		   "radius" : 1E-3,             // Radius		
		   "num particles": 1000,       // Number of computational particles
		   "weight" : 1                 // Particle weights
		} 	     
	     }
	  ]
       }
    ]

List/file
^^^^^^^^^

Particles can be read from a file by including a ``list`` specifier.
The particles must be organized as rows containing the coordinate positions and weights, e.g.

.. code-block::

   # x   y   z    w
     0   0   0    1

The user must specify

* ``file`` File name.
* ``x`` Column containing the :math:`x` coordinates (optional, defaults to 0).
* ``y`` Column containing the :math:`y` coordinates (optional, defaults to 1).
* ``z`` Column containing the :math:`z` coordinates (optional, defaults to 2).
* ``w`` Column containing the particle weight (optional, defaults to 3).

For example:

.. code-block:: json
   :emphasize-lines: 11-17
		     
    "plasma species" :
    [
       {
          "id": "e",                            // Species ID
          "Z" : 1,                              // Charge number
          "solver" : "ito",                     // Solver type. 
          "mobile" : false,                     // Not mobile
          "diffusive" : false,                  // Not diffusive
	  "initial particles": [                // Initial particles
	     {
	        "list": {
	   	   "file": "initial_particles.dat", // File containing the particles
		   "x column": 0,                   
		   "y column": 1,                   
		   "z column": 2,                   
		   "w column": 3
		}
	     }    
	  ]
       }
    ]

Initial densities
_________________

One may include an initial density by setting the ``initial density`` parameter.
For example:

.. code-block:: json
		     
    "plasma species" :
    [
       {
          "id": "e",           
          "Z" : 1,             
          "solver" : "cdr",    
          "mobile" : true,     
          "diffusive" : true,  
	  "initial density": 1E10
      }
    ]

If the species is defined as a particle species, computational particles will be generated within each grid so that the density is approximately as specified.

Photon species
--------------

Photon species are species that are tracked using a radiative transfer solver.
These are defined by an array of entries in a ``photon species`` JSON entry, and must define the following information:

* An ID/name for the species.
* The absorption coefficient for the photon type.

Basic definition
________________

A basic definition of a single photon species is

.. code-block:: json

   "photon species":
    [
	{
	    "id": "Y",              // Photon species id. Must be unique
	    "kappa": {              // Specification of absorption coefficient. 
		"type": "constant", // Specify constant absorption coefficient
		"value": 1E4        // Value of the absorption coefficeint
	    }
	}
    ]

Multiple photon species are added by appending with more entries, e.g.

.. code-block:: json

   "photon species":
    [
	{
	    "id": "Y1",             // Photon species id. Must be unique
	    "kappa": {              // Specification of absorption coefficient. 
		"type": "constant", // Specify constant absorption coefficient
		"value": 1E4        // Value of the absorption coefficeint
	    }
	},
	{
	    "id": "Y2",             // Photon species id. Must be unique
	    "kappa": {              // Specification of absorption coefficient. 
		"type": "constant", // Specify constant absorption coefficient
		"value": 1E4        // Value of the absorption coefficeint
	    }
	}	
    ]


Absorption coefficient
______________________

Constant
^^^^^^^^

In order to set a constant absorption coefficient, set ``type`` to constant and then specify the ``value`` field.
For example:

.. code-block:: json

   "photon species":
    [
	{
	    "id": "Y",
	    "kappa": { 
		"type": "constant",
		"value": 1E4       
	    }
	}
    ]

stochastic A
^^^^^^^^^^^^

The ``stochastic A`` specifier computes a random absorption length from the expression

.. math::

   \kappa = \left(\chi_{\textrm{min}}p_i\right)\left(\frac{\chi_{\textrm{max}}}{\chi_{\textrm{min}}}\right)^{\frac{f-f_1}{f_2 - f_1}},

where :math:`p_i` is the partial pressure of some species :math:`i`, and :math:`p_i\chi_{\textrm{min}}` and :math:`p_i\chi_{\textrm{max}}` are minimum and maximum absorption lengths on the frequency interval :math:`f\in[f_1,f_2]`.
The user must specify :math:`\chi_{\textrm{min}}`, :math:`\chi_{\textrm{max}}`, :math:`f_1`, :math:`f_2`, and the background species used when computing :math:`p_i`.
This is done through fields ``f1``, ``f_2``, ``chi min``, ``chi max``, and ``neutral``.
For example:

.. code-block:: json

    "photon species":
    [
	{
	    "id": "Y",                  
	    "kappa": {                  
		"type": "stochastic A", 
		"f1":   2.925E15,       
		"f2":   3.059E15,       
		"chi min": 2.625E-2,    
		"chi max": 1.5,         
		"neutral": "O2"         
	    }
	}
    ]		

stochastic B
^^^^^^^^^^^^

The ``stochastic B`` specifier computes a random absorption length from the expression

.. math::

   \kappa = \left(\chi_{\textrm{min}}p_i\right)\left(\frac{\chi_{\textrm{max}}}{\chi_{\textrm{min}}}\right)^u,

where :math:`p_i` is the partial pressure of some species :math:`i`, and :math:`p_i\chi_{\textrm{min}}` and :math:`p_i\chi_{\textrm{max}}` are minimum and maximum absorption lengths, and :math:`u` is a random number between 0 and 1.
Note that that this is just a simpler way of using the ``stochastic A`` specifier above. 
The user must specify :math:`\chi_{\textrm{min}}`, :math:`\chi_{\textrm{max}}`, and the background species used when computing :math:`p_i`.
This is done through fields ``chi min``, ``chi max``, and ``neutral``.
For example:

.. code-block:: json

    "photon species":
    [
	{
	    "id": "Y",                  
	    "kappa": {                  
		"type": "stochastic B", 
		"chi min": 2.625E-2,    
		"chi max": 1.5,         
		"neutral": "O2"         
	    }
	}
    ]		


Plasma reactions
----------------

Plasma reactions are always specified stoichiometrically in the form :math:`S_A + S_B + \ldots \xrightarrow{k} S_C + S_D + \ldots`.
The user must supply the following information:

* Species involved on the left-hand and right-hand side of the reaction.
* How to compute the reaction rate.
* Optionally specify gradient corrections to the reaction.
* Optionally specify whether or not the reaction rate will be written to file.

.. important::

   The JSON interface expects that the reaction rate :math:`k` is the same as in the reaction rate equation.
   Internally, the rate is converted to a format that is consistent with the KMC algorithm.

Plasma reactions should be entered in the JSON file using an array (the order of reactions does not matter).

.. raw:: html

   <details>
   <summary><a>Basic reaction example specifiers</a></summary>

.. code-block:: json

   "plasma reactions": [
      {
      // First reaction goes here,
      },
      {
      // Next reaction goes here,
      } 
   ]

.. raw:: html

   </details><br>      

Reaction specifiers
___________________

Each reaction must include an entry ``reaction`` in the list of reactions.
This entry must consist of a single string with a left-hand and right-hand sides separated by ``->``.
For example:

.. code-block:: json

   "plasma reactions": [
      {
         "reaction": "e + N2 -> e + e + N2+"
      }
   ]

which is equivalent to the reaction :math:`\text{e} + \text{N}_2 \rightarrow \text{e} + \text{e} + \text{N}_2^+`.
Each species must be separated by whitespace and a addition operator.
For example, the specification ``A + B + C`` will be interpreted as a reaction :math:`\text{A} + \text{B} + \text{C}` but the specification ``A+B + C`` will be interpreted as a reaction between one species ``A+B`` and another species ``C`` .

Internally, both the left- and right-hand sides are checked against background and plasma species.
The program will perform a run-time exit if the species are not defined in these lists.
Note that the right-hand side can additonal contain photon species.

.. tip::

   Products not not defined as a species can be enclosed by parenthesis ``(`` and ``)``.

The reaction string generally expects all species on each side of the reaction to be defined.
It is, however, occasionally useful to include products which are *not* defined.
For example, one may with to include the reaction :math:`\text{e} + \text{O}_2^+ \rightarrow \text{O} +\text{O}` but not actually track the species :math:`\text{O}`.
The reaction string can then be defined as the string ``e + O2+ ->`` but as it might be difficult for users to actually remember the full reaction, we may also enter it as ``e + O2+ -> (O) + (O)``. 

Rate calculation
________________

Reaction rates can be specified using multiple formats (constant, tabulated, function-based, etc).
To specify how the rate is computed, one must specify the ``type`` keyword in the JSON entry.

Constant
^^^^^^^^

To use a constant reaction rate, the ``type`` specifier must be set to constant and the reaction rate must be specified through the ``value`` keyword.
A JSON specification is e.g.

.. code-block:: json

   "plasma reactions": [
      {
         "reaction": "e + N2 -> e + e + N2+" // Specify reaction
	 "type": "constant",                 // Reaction rate is constant
	 "value": 1.E-30                     // Reaction rate
      }
   ]

Function vs E/N
^^^^^^^^^^^^^^^

Some rates given as a function :math:`k = k\left(E/N\right)` are supported, which are outlined below.

.. important::

   In the below function-based rates, :math:`E/N` indicates the electric field in units of Townsend.

**function E/N exp A**

This specification is equivalent to a fluid rate

.. math::

   k = c_1\exp\left[-\left(\frac{c_2}{c_3 + c_4 E/N}\right)^{c_5}\right].

The user must specify the constants :math:`c_i` in the JSON file.
An example specification is

.. code-block:: json
		
   "plasma reactions": [
      {
         "reaction": "A + B -> "       // Example reaction string. 
	 "type": "function E/N exp A", // Function based rate.
	 "c1": 1.0,                    // c1-coefficient
	 "c2": 1.0,                    // c2-coefficient
	 "c3": 1.0,                    // c3-coefficient
	 "c4": 1.0,                    // c4-coefficient
	 "c5": 1.0                     // c5-coefficient
      }
   ]


Temperature-dependent
^^^^^^^^^^^^^^^^^^^^^

Some rates can be given as functions :math:`k = k(T)` where :math:`T` is some temperature.

**function T A**

This specification is equivalent to a fluid rate

.. math::

   k = c_1\left(T_i\right)^{c_2}.

Mandatory input variables are :math:`c_1, c_2`, and the specification of the species corresponding to :math:`T_i`.
This can correspond to one of the background species.
An example specification is

.. code-block:: json
		
   "plasma reactions": [
      {
         "reaction": "A + B -> " // Example reaction string. 
	 "type": "function T A", // Function based rate.
	 "c1": 1.0,              // c1-coefficient
	 "c2": 1.0,              // c2-coefficient
	 "T": "A"               // Which species temperature
      }
   ]

**function T1T2 A**

This specification is equivalent to a fluid rate

.. math::

   k = c_1\left(\frac{T_1}{T_2}\right)^{c_2}.

Mandatory input variables are :math:`c_1, c_2`, and the specification of the species corresponding to :math:`T_1` and :math:`T_2`.
This can correspond to one of the background species.
An example specification is

.. code-block:: json
		
   "plasma reactions": [
      {
         "reaction": "A + B -> "    // Example reaction string. 
	 "type": "function T1T2 A", // Function based rate.
	 "c1": 1.0,                 // c1-coefficient
	 "c2": 1.0,                 // c2-coefficient
	 "T1": "A",                 // Which species temperature for T1
	 "T2": "B"                  // Which species temperature for T2	 
      }
   ]   

Townsend rates
^^^^^^^^^^^^^^

Reaction rates can be specified to be proportional to :math:`\alpha\left|\mathbf{v}_i\right|` where :math:`\alpha` is the Townsend ionization coefficient and :math:`\left|\mathbf{v}_i\right|` is the drift velocity for some species :math:`i`.
This type of reaction is normally encountered when using simplified chemistry, e.g.

.. math::

   \partial_t n_i = \alpha\left|\mathbf{v}_s\right| n_i = \alpha\mu\left|E\right|n_i.

which is representative of the reaction :math:`S_i \rightarrow S_i + S_i`.
To specify a Townsend rate constant, one can use the following:

#. ``alpha*v`` for setting the rate constant proportional to the Townsend ionization rate.
#. ``eta*v`` for setting the rate constant proportional to the Townsend attachment rate.

One must also specify which species is associated with :math:`\left|\mathbf{v}\right|` by specifying a species flag.
A complete JSON specification is

..
   .. raw:: html

      <details>
      <summary><a>Example Townsend rate JSON specification</a></summary>

.. code-block:: json

   {
      "plasma reactions":
         [
	    {
	       "reaction": "e -> e + e + M+", // Reaction string
	       "type": "alpha*v",             // Rate is alpha*v
	       "species": "e"                 // Species for v
	    }	
         ]
   }

..
   .. raw:: html

      </details><br>   

.. warning::

   When using the Townsend coefficients for computing the rates, one should normally *not* include any neutrals on the left hand side of the reaction.
   The reason for this is that the Townsend coefficients :math:`\alpha` and :math:`\eta` already incorporate the neutral density.
   By specifying e.g. a reaction string :math:`\text{e} + \text{N}_2 \rightarrow \text{e} + \text{e} + \text{N}_2^+` together with the ``alpha*v`` or ``eta*v`` specifiers, one will end up multiplying in the neutral density twice, which will lead to an incorrect rate.

Tabulated vs E/N
^^^^^^^^^^^^^^^^   

Scaling
_______

Reactions can be scaled by including a ``scale`` field in the JSON entry.
This will scale the reaction coefficient by the input factor, e.g. modify the rate constant as

.. math::

   k \rightarrow \nu k,
   
where :math:`\nu` is the scaling factor.
This is useful when scaling reactions from different units, or for completely turning off some input reactions.
An example JSON specification is

.. code-block:: json
		
    "plasma reactions":
    [
	{
	    "reaction": "e -> e + e + M+", // Reaction string
	    "type": "alpha*v",             // Rate is alpha*v
	    "species": "e",                // Species for v,
	    "scale": 0.0                   // Scaling factor
	}	
    ]

.. tip::

   If one turns off a reaction by setting ``scale`` to zero, the KMC algorithm will still use the reaction but no reactants/products are consumed/produced.

Efficiency
__________

Reactions efficiencies can be modified in the same way as one do with the ``scale`` field, e.g. modify the rate constant as

.. math::

   k \rightarrow \nu k,

where :math:`\nu` is the reaction efficiency.
Specifications of the efficiency can be achieved in the forms discussed below.

Constant efficiency
^^^^^^^^^^^^^^^^^^^

An example JSON specification for a constant efficiency is:

.. code-block:: json
		
    "plasma reactions":
    [
	{
	    "reaction": "e -> e + e + M+", 
	    "efficiency": 0.5              
	}	
    ]

Efficiency vs E/N
^^^^^^^^^^^^^^^^^

An example JSON specification for an efficiency computed versus :math:`E/N` is

.. code-block:: json
		
    "plasma reactions":
    [
	{
	    "reaction": "e -> e + e + M+", 
	    "efficiency vs E/N": "efficiencies.dat"
	}	
    ]

where the file ``efficiencies.dat`` must contain two-column data containing values of :math:`E/N` along the first column and efficiencies along the second column.
An example file is e.g.

.. code-block:: text

   0    0
   100  0.1
   200  0.3
   500  0.8
   1000 1.0
   
This data is then internally convered to a uniformly spaced lookup table (see :ref:`Chap:LookupTable`).

Efficiency vs E
^^^^^^^^^^^^^^^^^

An example JSON specification for an efficiency computed versus :math:`E` is

.. code-block:: json
		
    "plasma reactions":
    [
	{
	    "reaction": "e -> e + e + M+", 
	    "efficiency vs E": "efficiencies.dat"
	}	
    ]

where the file ``efficiencies.dat`` must contain two-column data containing values of :math:`E` along the first column and efficiencies along the second column.
This method follows the same as ``efficiency vs E/N`` where the data in the input file is put in a lookup table.

Quenching
_________

Reaction quenching can be achieved in the following forms:

Pressure-based
^^^^^^^^^^^^^^

The reaction rate can modified by a factor :math:`p_q/\left[p_q + p\left(\mathbf{x}\right)\right]` where :math:`p_q` is a quenching pressure and :math:`p\left(\mathbf{x}\right)` is the gas pressure.
This will modify the rate as

.. math::

   k \rightarrow k\frac{p_q}{p_q + p\left(\mathbf{x}\right)}.

An example JSON specification is

.. code-block:: json
		
    "plasma reactions":
    [
	{
	    "reaction": "e + N2 -> e + N2 + Y", // Reaction string
	    "type": "alpha*v",                  // Rate is alpha*v
	    "species": "e",                     // Species for v,
	    "quenching pressure": 4000.0        // Quenching pressure
	}	
    ]

.. important::

   The quenching pressure must be specified in units of Pascal.

Rate-based
^^^^^^^^^^

The reaction rate can be modified by a factor :math:`k_r/\left[k_r + k_p + k_q\right]`.
The intention behind this scaling is that reaction :math:`r` occurs only if it is not predissociated (by rate :math:`k_p`) or quenched (by rate :math:`k_q`).
Such processes can occur, for example, in excited molecules.
This will modifiy the rate constant as

.. math::

   k \rightarrow k\frac{k_r}{k_r + k_p + k_q},

where the user will specifiy :math:`k_r`, :math:`k_p`, and :math:`k_q/N`.
The JSON specification must contain ``quenching rates``, for example:

.. code-block:: json
		
    "plasma reactions":
    [
	{
	    "reaction": "e + N2 -> e + N2 + Y", // Reaction string
	    "type": "alpha*v",                  // Rate is alpha*v
	    "species": "e",                     // Species for v,
	    "quenching rates": {                // Specify relevant rates
	       "kr": 1E9,
	       "kp": 1E9,
	       "kq/N": 0.0 
	    }
	}	
    ]




Gradient correction
___________________

In LFA-based models it is frequently convenient to include energy-corrections to the ionization rate.
In this case one modifies the rate as

.. math::

   k \rightarrow k\left(1 + \frac{\mathbf{E}\cdot \left(D\nabla n_i\right)}{n_i\mu_iE^2}\right).

To include this correction one may include a specifier ``gradient correction`` in the JSON entry, in which case one must also enter the species :math:`n_i`.
A JSON specification that includes this

.. code-block:: json
		
    "plasma reactions":
    [
	{
	    "reaction": "e -> e + e + M+", // Reaction string
	    "type": "alpha*v",             // Rate is alpha*v
	    "species": "e",                // Species for v,
	    "gradient correction": "e"     // Specify gradient correction using species "e"
	}	
    ]

This is equivalent to a source term

.. math::

   S &= \alpha n_e\left|\mathbf{v}_e\right|\left(1 + \frac{\mathbf{E}\cdot \left(D_e\nabla n_e\right)}{n_e\mu_eE^2}\right) \\
   &=\alpha\left[n_e\left|\mathbf{v}_e\right| + \hat{\mathbf{E}}\cdot\left(D_e\nabla n_e\right)\right].

One can recognize this term as a regular electron impact ionization source term (typically written as :math:`\alpha \mu n_e E`).
With the gradient correction, the ionization source term is essentially computed using the full electron flux, i.e., including the diffusive electron flux.
Note that the full electron flux has a preferential direction, and the physical interpretation of this direction is that if there is net diffusion against the electric field, electrons lose energy and the impact ionization source term is correspondingly lower.



Understanding reaction rates
____________________________

The JSON specification takes *fluid rates* as input to the KMC algorithm.
Note that subsequent application of multiple scaling factors are multiplicative.
For example, the specification

.. code-block:: json
		
    "plasma reactions":
    [
	{
	    "reaction": "e -> e + e + M+", // Reaction string
	    "type": "alpha*v",             // Rate is alpha*v
	    "species": "e",                // Species for v,
	    "gradient correction": "e"     // Specify gradient correction using species "e"
	    "quenching pressure": 4000,    // Quenching pressure
	    "efficiency": 0.1              // Reaction efficiency
	}	
    ]

is equivalent to the rate constant

.. math::

   k = \alpha\left|\mathbf{v}_e\right|\left(1 + \frac{\mathbf{E}\cdot \left(D_e\nabla n_e\right)}{n_e\mu_eE^2}\right)\frac{p_q}{p_q+p\left(\mathbf{x}\right)}\nu

Internally, conversions between the *fluid rate* :math:`k` and the KMC rate :math:`c` are made.
The conversion factors depend on the reaction order, and ensure consistency between the KMC and fluid formulations.

Plotting rates
______________

Reaction rates can be plotted by including an optional ``plot`` specifier.
If the specifier included, the reaction rates will be included in the HDF5 output.
The user can also include a description string which will be used when plotting the reaction.

The following two specifiers can be included:

* ``plot`` (true/false) for plotting the reaction.
* ``description`` (string) for naming the reaction in the HDF5 output.

If the plot description is left out, the reaction string will be used as a variable name in the plot file.
A JSON specification that includes these

.. code-block:: json
		
    "plasma reactions":
    [
	{
	    "reaction": "e -> e + e + M+",       // Reaction string
	    "type": "alpha*v",                   // Rate is alpha*v
	    "species": "e",                      // Species for v,
	    "plot": true,                        // Plot this reaction
	    "description": "Townsend ionization" // Variable name in HDF5
	}	
    ]

Printing rates
______________

``ItoKMCJSON`` can print the plasma reaction rates (including photon-generating ones) to file.
This is done through an optional input argument ``print_rates`` in the input file (not the JSON specification).
The following options are supported:

.. code-block:: text

   ItoKMCJSON.print_rates            = true
   ItoKMCJSON.print_rates_minEN      = 1
   ItoKMCJSON.print_rates_maxEN      = 1000
   ItoKMCJSON.print_rates_num_points = 200
   ItoKMCJSON.print_rates_spacing    = exponential
   ItoKMCJSON.print_rates_filename   = fluid_rates.dat
   ItoKMCJSON.print_rates_pos        = 0 0 0

Setting ``ItoKMCJSON.print_rates`` to true in the input file will write all reaction rates as column data :math:`E/N, k(E/N)`.
Here, :math:`k` indicates the *fluid rate*, so for a reaction :math:`A + B + C \xrightarrow{k}\ldots` it will include the rate :math:`k`.
Reactions are ordered identical to the order of the reactions in the JSON specification.
This feature is mostly used for debugging or development efforts.

Particle placement
------------------

Specification of secondary particle placement is done through the JSON file by specifying the ``particle placement`` field.
Currently, new particles may be placed on the centroid, uniformly distributed in the grid cell, or placed randomly in the downstream region of some user-defined species.
The method is specified through the ``method`` specifier, which must either be ``centroid``, ``random``, or ``downstream``.
If specifying ``downstream``, one must also include a species specifier (which must correspond to one of the plasma species).
The following three specifies are all valid:

.. code-block:: json
		
    "particle placement":
    {
       "method": "centroid"
    }

    "particle placement":
    {
       "method": "random"
    }

    "particle placement":
    {
       "method": "downstream",
       "species": "e"
    }        

Photoionization
---------------

Photoionization reactions are specified by including a ``photoionization`` array of JSON entries that specify each reaction.
The left-hand side of the reaction must contain a photon species (plasma and background species can be present but will be ignored), and the right-hand side can consist of any species except other photon species.
Each entry *must* contain a reaction string specifier and optionally also an efficiency (which defaults to 1).
Assume that photons are generated through the reaction

.. code-block:: json
		
    "plasma reactions":
    [
	{
	    "reaction": "e + N2 -> e + N2 + Y", // Reaction string
	    "type": "alpha*v",                  // Rate is alpha*v
	    "species": "e",                     // Species for v,
	    "quenching pressure": 4000.0        // Quenching
	    "efficiency": 0.6                   // Excitation events per ionization event
	    "scale": 0.1                        // Photoionization events per absorbed photon
	}	
    ]

where we assume an excitation efficiency of :math:`0.6` and photoionization efficiency of :math:`0.1`.
An example photoionization specification is then

.. code-block:: json
		
    "photoionization":
    [
	{
	    "reaction": "Y + O2 -> e + O2+"
	}	
    ]


Note that the efficiency of this reaction is 1, i.e. the photoionization probabilities was pre-evaluated.
As discussed in :ref:`Chap:ItoKMCPhysics`, we can perform late evaluation of the photoionization probability by specifically including a efficiency.
In this case we modify the above into

.. code-block:: json

    "plasma reactions":
    [
	{
	    "reaction": "e + N2 -> e + N2 + Y", // Reaction string
	    "type": "alpha*v",                  // Rate is alpha*v
	    "species": "e",                     // Species for v,
	    "quenching pressure": 4000.0        // Quenching
	    "efficiency": 0.6                   // Excitation events per ionization event
	}	
    ],
    "photoionization":
    [
	{
	    "reaction": "Y + O2 -> e + O2+",
	    "efficiency": 0.1
	},
	{
	    "reaction": "Y + O2 -> (null)",
	    "efficiency": 0.9
	}		
    ]

where a null-absorption model has been added for the photon absorption.
When multiple pathways are specified this way, and they have probabilities :math:`\xi_1, \xi_2, \ldots`, the reaction is stochastically determined from a discrete distribution with relative probabilities :math:`p_i = \xi_i/(\xi_1 + \xi_2+\ldots)`.

.. important::

   The null-reaction model is *not* automatically added when using late evaluation of the photoionization probability.


Secondary emission
------------------

Secondary emission is supported for particles, including photons, but not for fluid species.
Specification of secondary emission is done separately on dielectrics and electrodes by specifying JSON entries ``electrode emission`` and ``dielectric emission``.
The internal specification for these are identical, so all examples below use ``dielectric emission``.
Secondary emission reactions are specified similarly to photoionization reactions.
However, the left-hand side must consist of a single species (photon or particle), while the right-hand side can consist of multiple species.
Wildcards ``@`` are supported, as is the ``(null)`` species which enables null-reactions, as further discussed below.
An example JSON specification that enables secondary electron emission due to photons and ion impact is

.. code-block:: json

   "dielectric emission":
   [
      {
         "reaction": "Y -> e",
	 "efficiency": 0.1
      },
      {
         "reaction": "Y -> (null)",
	 "efficiency": 0.9
      }
   ]

The ``reaction`` field specifies which reaction is triggered: The left hand side is the primary (outgoing) species and the left hand side must contain the secondary emissions.
The left-hand side can consist of either photons (e.g., ``Y``) or particle (e.g., ``O2+``) species.
The efficiency field specifies the efficiency of the reaction.
When multiple reactions are specified, we randomly sample the reaction according to a discrete distribution with probabilities

.. math::

   p(i) = \frac{\nu_i}{\sum_i \nu_i},

where :math:`\nu_i` are the efficiencies.
Note that the efficiencies do not need to sum to one, and if only a single reaction is specified the efficiency specifier has not real effect.
The above reactions include the null reaction in order to ensure that the correct secondary emission probability is reached, where the ``(null)`` specifier implies that no secondary emission takes place.
In the above example the probability of secondary electron emission is 0.1, while the probability of a null-reaction (outgoing particle is absorbed without any associated emission) is 0.9.
The above example can be compressed by using a wildcard and an ``efficiencies`` array as follows:

.. code-block:: json

   "dielectric emission":
   [
      {
      "reaction": "Y -> @",
	 "@": ["e", "(null)"],
	 "efficiencies": [1,9]
      }
   ]

where for the sake of demonstration the efficiencies are set to 1 and 9 (rather than 0.1 and 0.9).
This has no effect on the probabilities :math:`p(i)` given above.

Field emission
--------------

Field emission for a specified species is supported by including a ``field emission`` specifier.
Currently supported expressions are:

#. Fowler-Nordheim emission:

   .. math::

      J(E) &= \frac{a F^2}{\phi}\exp\left(-v(f) b \phi^{3/2} F\right),\\
      F &= \left(\beta E\right)\times 10^9,\\
      a &= 1.541434\times 10^{-6}\,\text{AeV}/\text{V}^2, \\
      b &= 6.830890\,\text{eV}^{-3/2}\text{V}\text{nm}^{-1},\\
      v(f) &= 1 - f + \left(1/6\right)f\log(f),\\
      f &\approx 1.439964\,\text{eV}^2\text{V}^{-1}\text{nm}\times \frac{F}{\phi^2}.

   Here, :math:`\phi` is the work function (in electron volts) and :math:`\beta` is an empirical factor that describes local field amplifications.
   Note that the above expression gives :math:`J` in units of :math:`\text{A}/\text{nm}^2`.
      

#. Schottky emission:

   .. math::

      J(E) &= \lambda A_0 T^2\exp\left(-\frac{\left(\phi-\Delta\phi\right)q_\text{e}}{k_{\text{B}}T}\right), \\
      F &= \beta E, \\
      A_0 &\approx 1.201736\times 10^6\,\text{A}\text{m}^{-2}\text{K}^{-2}, \\
      \Delta \phi &= \sqrt{\frac{q_\text{e} F}{4\pi\epsilon_0}}.

   Here, :math:`T` is the cathode temperature and :math:`\phi` is the work function (in electron volts).
   The value :math:`\lambda` is typically around :math:`0.5`.
   As for the Fowler-Nordheim equation, a factor that describes local field amplifications is given in :math:`\beta`. 

.. warning::
      
   The expressions above both use a correction factor :math:`\beta` that describes local field amplifications.
   Note, however, that the number of emitted electrons is proportional to the area of the emitter, and there are no current models in ``chombo-discharge`` that can compute this effective area.

Below, we show an example that includes both Schottky and Fowler-Nordheim emission

.. code-block:: json

    "field emission":
    [
	{
	    "species": "e",
	    "surface": "electrode",
	    "type": "fowler-nordheim",
	    "work": 5.0,
	    "beta": 1
	},
	{
	    "species": "e",
	    "surface": "electrode",
	    "type": "schottky",
	    "lambda": 0.5	    
	    "temperature": 300,
	    "work": 5.0,
	    "beta": 1,
	}
    ]		


.. _Chap:ItoKMCWarnings:
		
Warnings and caveats
--------------------

Higher-order reactions
______________________

Usually, many rate coefficients depend on the output of other software (e.g., BOLSIG+) and the scaling of rate coefficients is not immediately obvious.
This is particularly the case for three-body reactions with BOLSIG+ that may require scaling before running the Boltzmann solver (by scaling the input cross sections), or after running the Boltzmann solver, in which case the rate coefficients themselves might require scaling. In any case the user should investigate the cross-section file that BOLSIG+ uses, and figure out the required scaling. 

.. important::
   
   For two-body reactions, e.g. :math:`A + B\rightarrow \varnothing` the rate coefficient must be specified in units of :math:`\textrm{m}^3\textrm{s}^{-1}`, while for three-body reactions :math:`A + B + C\rightarrow\varnothing` the rate coefficient must have units of :math:`\textrm{m}^6\textrm{s}^{-1}`.

   For three-body reactions the units given by BOLSIG+ in the output file may or may not be incorrect (depending on whether or not the user scaled the cross sections).

Townsend coefficients
_____________________

Townsend coefficients are not fundamentally required for specifying the reactions, but as with the higher-order reactions some of the output rates for three-body reactions might be inconsistently represented in the BOLSIG+ output files.
For example, some care might be required when using the Townsend attachment coefficient for air when the reaction :math:`\textrm{e} + \textrm{O}_2 + \textrm{O}_2\rightarrow \textrm{O}_2^- + \textrm{O}_2` is included because the rate constant might require proper scaling after running the Boltzmann solver, but this scaling is invisible to the BOLSIG+'s calculation of the attachment coefficient :math:`\eta/N`.

.. warning::

   The JSON interface *does not guard* against inconsitencies in the user-provided chemistry, and provision of inconsistent :math:`\eta/N` and attachment reaction rates are quite possible.

Tips and tricks
---------------

As with fluid drift-diffusion models, numerical instabilities can also occur due to unbounded growth in the plasma density.
This is a process which has been linked both to the local field approximation and also to the presence of numerical diffusion.
Simulations that fail to stabilize, i.e., where the field strength diverges, may benefit from the following stabilizing features:

#. **Turn off parallel diffusion.**
   
   With the ``ItoKMCGodunovStepper`` class, this option is given by ``ItoKMCGodunovStepper.limit_parallel_diffusion``.
   Using this option will ensure that particles do not diffuse against their drift direction.
   Note that this also modifies the amount of *physical* diffusion in this direction.

#. **Use gradient corrections.**

   As discussed earlier, using a gradient correction can help limit non-physical ionization due to backwards-diffusing electrons.

#. **Use downstream particle placement.**

   Because the KMC algorithm solves for the number of particles in a grid cell, distributing new particles uniformly over a grid cell can lead to numerical diffusion where secondary electrons are placed in the wake of primary electrons.
   Using downstream particle placement often leads to slightly more stable simulations.

#. **Describe the primary species using an Ito solver.**

   Similar to the point above, using a fluid solver for the ions may lead to upstream placement of the resulting positive charge.   

#. **Use kd-tree particle merging.**

   Currently particle merging strategies are reinitialization and merging based on bounding volume hierarchies.
   If using reinitialization, new particles can be generated in the wake of old ones and can thus upset the charge distribution and cause numerical backwards diffusion (e.g., of electrons).

#. **Numerically limit reaction rates.**

   It is possible to specify that reaction rates will be numerically limited so that the rate does not exceed a specified threshold of :math:`\Delta t^{-1}`.
   This is done through the JSON interface by setting ``limit max k*dt`` to some value.
   Note that this changes the physics of the model, but usually enhances stability at larger time steps.
   An example is given below:

.. code-block:: json

   "plasma reactions": [
      {
         "reaction": "e + N2 -> e + e + N2+"
	 "type": "constant",
	 "value": 1.E-12,
	 "limit max k*dt" : 2.0
      }
   ]

This will limit the rate such that :math:`k \left[\text[N]_2\right]\Delta t = 2`.
I.e., all background species are first absorbed into the rate calculation before the rate is limited.
We point out that limiting is not possible if both species on the left hand side are solver variables.


.. important::
   
   The above features have been implemented in order to push the algorithm towards coarser grids and larger time steps.
   It is essential that the user checks that the model converges when these features are applied.
      
Example programs
================

Example programs that use the Îto-KMC module are given in

* :file:`$DISCHARGE_HOME/Exec/Examples/ItoKMC/AirBasic` for a basic streamer discharge in atmospheric air.
* :file:`$DISCHARGE_HOME/Exec/Examples/ItoKMC/AirDBD` for a streamer discharge over a dielectric.
