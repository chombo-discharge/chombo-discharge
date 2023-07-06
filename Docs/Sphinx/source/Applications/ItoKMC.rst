.. _Chap:KMC:

Îto-KMC plasma model
********************

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

Reactions
---------

Particle management
-------------------

0D chemistry
------------

The user input interface to the Îto-KMC model consists of a zero-dimensional plasma kinetics interface called ``ItoKMCPhysics``.
This interface consists of the following main functionalities:

#. User-defined species, and which solver types (Îto or CDR) to use when tracking them.
#. User-defined initial particles, reaction kinetics, and photoionization processes. 

The complete interface specification is given below.
Because the interface is fairly extensive, ``chombo-discharge`` also supplies a JSON-based implementation called ``ItoKMCJSON`` (see :ref:`Chap:ItoKMCJSON`) for defining these things through file input.
The difference between ``ItoKMCPhysics`` and its implementation ``ItoKMCJSON`` is that the JSON-based implementation class only supports a subset of potential features supported by ``ItoKMCPhysics``.

.. raw:: html

   <details>
   <summary><a>Show/hide full 0D interface</a></summary>

.. literalinclude:: ../../../../Physics/ItoKMC/CD_ItoKMCPhysics.H
   :language: c++   

.. raw:: html

   </details><br>

.. _Chap:ItoKMCJSON:

JSON 0D chemistry interface
===========================

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

Each user-defined gas law exists as a separate entry in the ``gas/law`` field, and the minimum requirements for these are that they contain a ``type`` specifier.
Currently, the supported types are

* ``ideal`` (for constant ideal gas)
* ``table vs height`` for tabulated density, temperature and pressure vs some height/axis.
  
The specification rules for these are different, and are discussed below.
An example JSON specification is given below.
Here, we have defined two gas law ``my_ideal_gas`` and ``tabulated_atmosphere`` and specified through the ``id`` parameter that we will use the ``my_ideal_gas`` specification.

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
		// Definition for an ideal gas law. Temperature is always in Kevlin the pressure is in bar. The neutral density
		// is derived using an ideal gas law.
		"type" : "ideal",
		"temperature" : 300,
		"pressure" : 1.0
	    },
	    "tabulated_atmosphere" : {
		// Tabulated gas law. The user must supply a file containing the pressure, temperature and number density and
		// specify the table resolution that will be used internally.
		"type" : "table vs height",       // Specify that our gas law contains table-vs-height data
		"file" : "ENMSIS_Atmosphere.dat", // File name
		"axis" : "y",                     // Associated Cartesian axis for the "height"
		"temperature column" : 0,         // Column containing the temperature
		"pressure column" : 1,            // Column containing the pressure
		"density column" : 2,             // Column containing the number density
		"min height" : 0.0,               // Minium height kept when resampling the table
		"max height" : 250000,            // Maximum height kept when resampling the table
		"res height" : 500,               // Table resolution
		"dump tables" : true,             // Optional argument for dumping tables to disk (useful for debugging)
		"T scale" : 1.0,                  // Optional argument for user-specified temperature data scaling.
		"P scale" : 1.0,                  // Optional argument for user-specified pressure data scaling.
		"Rho scale" : 1.0                 // Optional argument for user-specified density data scaling.
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

.. _Chap:ItoKMCJSON:

Table vs height
_______________


Background species
------------------

Name
____

Molar fraction
______________

Plotting
________

Townsend coefficients
---------------------

Plasma species
--------------

Basic definition
________________

Initial particles
_________________

Mobility coefficients
_____________________


Diffusion coefficients
______________________

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

Plasma reactions should be entered in the JSON file using an array, .e.g.


.. code-block:: json

   "plasma reactions": [
      {
      // First reaction goes here,
      },
      {
      // Next reaction goes here,
      } 
   ]

The order of reactions does not matter.

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

which is equivalent to the reaction :math:`\text{e} + \text{N}_2 \xrightarrow{k} \text{e} + \text{e} + \text{N}_2^+`.
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

**function TT A**

This specification is equivalent to a fluid rate

.. math::

   k = c_1\left(\frac{T_1}{T_2}\right)^{c_2}.

Mandatory input variables are :math:`c_1, c_2`, and the specification of the species corresponding to :math:`T_1` and :math:`T_2`.
This can correspond to one of the background species.
An example specification is

.. code-block:: json
		
   "plasma reactions": [
      {
         "reaction": "A + B -> "  // Example reaction string. 
	 "type": "function TT A", // Function based rate.
	 "c1": 1.0,               // c1-coefficient
	 "c2": 1.0,               // c2-coefficient
	 "T1": "A",               // Which species temperature for T1
	 "T2": "B"                // Which species temperature for T2	 
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
.. A complete JSON specification is

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

.. raw:: html

   </details><br>   

.. warning::

   When using the Townsend coefficients for computing the rates, one should normally *not* include any neutrals on the left hand side of the reaction.
   The reason for this is that the Townsend coefficients :math:`\alpha` and :math:`\eta` already incorporate the neutral density.
   By specifying e.g. a reaction string :math:`\text{e} + \text{N}_2 \rightarrow \text{e} + \text{e} + \text{N}_2^+` together with the ``alpha*v`` or ``eta*v`` specifiers, one will end up multiplying in the neutral density twice, which will lead to an incorrect rate.

Tabulated vs E/N
^^^^^^^^^^^^^^^^   

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

Scaling
_______

Reactions can be scaled by including a ``scale`` field in the JSON entry.
This will scale the reaction coefficient by the input factor.
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
      
Example programs
================
