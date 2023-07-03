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

JSON interface
==============

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

Rate calculation
________________

Constant
^^^^^^^^

Tabulated vs E/N
^^^^^^^^^^^^^^^^

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
   By specifying e.g. a reaction string :math:`\text{e} + N_2 \rightarrow \text{e} + \text{e} + N_2^+` together with the ``alpha*v`` or ``eta*v`` specifiers, one will end up multiplying in the neutral density twice, which will lead to an incorrect rate.







