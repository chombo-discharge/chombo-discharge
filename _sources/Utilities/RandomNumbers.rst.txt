.. _Chap:Random:

Random numbers
==============

``Random`` is a static class for generating pseudo-random numbers, and exist so that all random number operations can be aggregated into a single class.
Internally, ``Random`` use a Mersenne-Twister random number generation supplied by ``std::random``. 

To use the ``Random`` class, simply include ``<CD_Random.H>``, e.g.,

.. code-block:: c++

   #include <CD_Random.H>

See the `Random API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classRandom.html>`_ for further details. 

Drawing random numbers
----------------------

General approach
________________

The general routine for drawing a random number is

.. literalinclude:: ../../../../Source/Utilities/CD_Random.H
   :language: c++
   :lines: 111-117
   :dedent: 2

Here, the template parameter ``T`` is some distribution that follows the appropriate C++ template constraints of ``<random>``.
For example:

.. code-block:: c++

   std::uniform_real_distribution<Real> dist(0.0, 100.0);

   const Real randomNumber = Random::get(dist);

Pre-defined distributions
_________________________

Pre-defined distributions exist for performing the following operations:

.. literalinclude:: ../../../../Source/Utilities/CD_Random.H
   :language: c++
   :lines: 69-109
   :dedent: 2

Seeding the RNG
---------------

By default, the random number generator is seeded using the MPI rank.
It is necessary to seed MPI ranks using different seeds to avoid them producing the same number sequences.
If using MPI+OpenMP, additional steps are taken to ensure that each OpenMP thread obtains a unique random number generator.

``Driver`` (see :ref:`Chap:Driver`) will seed the random number generator, and user can override the seed by setting ``Random.seed = <number>`` in the input script.
If the user sets ``<number> < 0`` then a random seed will be produced based on the elapsed CPU clock time.
If running with MPI, this seed is obtained by only one of the MPI ranks, and this seed is then broadcast to all the other ranks.
The other ranks will then increment the seed by their own MPI rank number so that each MPI rank gets a unique seed.

