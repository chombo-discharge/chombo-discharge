.. _Chap:KMC:

Îto-KMC plasma model
********************

.. warning::

   Documentation is pending.

Underlying model
================


Algorithms
==========

Time stepping
-------------

Reactions
---------

Particle management
-------------------

JSON 0D chemistry interface
===========================

Gas law
-------

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
A complete JSON specification is

.. code-block:: json
		
    "plasma reactions":
    [
	{
	    "reaction": "e -> e + e + M+", // Reaction string
	    "type": "alpha*v",             // Rate is alpha*v
	    "species": "e"                 // Species for v
	}	
    ]

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
This is useful when scaling reactions from different units, or for completely off some input reactions.
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
