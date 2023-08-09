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

Reaction network
----------------

Particle management
-------------------

0D chemistry
------------

The user input interface to the Îto-KMC model consists of a zero-dimensional plasma kinetics interface called ``ItoKMCPhysics``.
This interface consists of the following main functionalities:

#. User-defined species, and which solver types (Îto or CDR) to use when tracking them.
#. User-defined initial particles, reaction kinetics, and photoionization processes.

.. _Chap:ItoKMCPhysics:

ItoKMCPhysics
_____________

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
___________________

Species are defined either as input species for CDR solvers (see :ref:`Chap:CdrSolver`) or Îto solvers (see :ref:`Chap:ItoSolver`).
It is sufficient to populate the ``ItoKMCPhysics`` species vectors ``m_cdrSpecies`` and ``m_itoSpecies`` if one only wants to initialize the solvers.
See :ref:`Chap:ItoKMCPhysics` for their C++ definition.
In addition to actually populating the vectors, users will typically also include initial conditions for the species.
The interfaces permit initial particles for the Îto solvers, while the CDR solvers permit particles *and* initial density functions.
Additionally, one must define all species associated with radiative transfer solvers.

.. _Chap:ItoKMCPlasmaReaction:

Plasma reactions
________________

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
Internally, these reactions are implemented through the dual state KMC implementation, see :ref:`Chap:KMCDualState`.
During the reaction advance the user only needs to update the :math:`c_j` coefficients (typically done via an interface implementation); the calculation of the propensity is automatic and follows the standard KMC rules (e.g., the KMC solver accounts for the number of distinct pairs of particles).
This must be done in the routine ``updateReactionRates(...)``, see :ref:`Chap:ItoKMCPhysics` for the complete specification.

Photo-reactions
_______________

Photo-reactions are also represented stoichiometrically as

.. math::

   \gamma \xrightarrow{\xi} A + B + \ldots,

where :math:`\xi` is the photo-reaction efficiency, i.e. the probability that the photon :math:`\gamma` causes a reaction when it is absorbed on the mesh.
Currently, photons can only be generated through *plasma reactions*, see :ref:`Chap:ItoKMCPlasmaReaction`, and we do not permit reactions between the photon and plasma species, i.e. reactions of the type :math:`\gamma + A\rightarrow\varnothing` are not permitted. 

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
Pre-scaling by the photo-reaction efficiency is then difficult because the reactions

.. math::

   \gamma &\xrightarrow{\text{volume}}\varnothing \\
   \gamma &\xrightarrow{\text{surface}}\varnothing,

can not be pre-scaled.

.. note::

   Is this true? Or can we scale by :math:`\sum_j\xi_j + \sum_j\zeta_j` and then evaluate :math:`\xi_i/\sum_j \xi_j` when doing absorption on in the volume?

   Should this be automated? How should the photon-products be scaled???
   
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

.. important::
   
   The above sampling routines require some modification when dealing with super-photons.
   For example, if the photon weight is :math:`w_\gamma = 1000` and :math:`\xi_1 = 0.01`, absorption of the super-photon on the mesh should on average trigger 10 reactions.
   

Surface reactions
_________________

Transport coefficients
______________________

Species mobilities and diffusion coefficients should be computed just as they are done in the fluid approximation.
The ``ItoKMCPhysics`` interface requires implementations of two functions that define the coefficients as functions of :math:`\mu = \mu\left(t,\mathbf{x}, \mathbf{E}\right)`,
and likewise for the diffusion coefficients.
Note that these functions should return the *fluid coefficients*.

.. important::

   There is currently no support for computing :math:`\mu` as a function of the species densities (e.g., the electron density), but this only requires modest extensions of the Îto-KMC module.

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


Table vs height
_______________

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

.. code-block:: txt

   # height       molar fraction
   0              0.1
   1              0.1 
   2              0.1

The file parser (see :ref:`LookupTable`) will ignore the header file if it starts with a hashtag (#).
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
   :emphasize-lines: 6-18
      
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
#. Usage in plasma reactions.

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

  .. code-block:: txt

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
An example JSON specification is

.. code-block:: json
		
    "plasma reactions":
    [
	{
	    "reaction": "e -> e + e + M+", // Reaction string
	    "type": "alpha*v",             // Rate is alpha*v
	    "species": "e",                // Species for v,
	    "efficiency": 0.5              // Reaction efficiency.
	}	
    ]

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

Photoionization
---------------
      
Example programs
================
