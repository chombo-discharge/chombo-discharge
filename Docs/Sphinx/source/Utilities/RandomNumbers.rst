.. _Chap:Random:

Random numbers
==============

``Random`` is a static class for generating pseudo-random numbers, and exist so that all random number operations can be aggregated into a single class.
Internally, ``Random`` use a Mersenne-Twister random number generation.

To use the ``Random`` class, simply include ``<CD_Random.H>``, e.g.

.. code-block:: c++

   #include <CD_Random.H>

See the `Random API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classRandom.html>`_ for further details. 

Drawing random numbers
----------------------

The general routine for drawing a random number is

.. code-block:: c++

  template<typename T>
  Real get(T& a_distribution);

which is for example used as follows:

.. code-block:: c++

   std::uniform_real_distribution<Real> dist(0.0, 100.0);

   const Real randomNumber = Random::get(dist);

Pre-defined distributions exist for performing the following operations:

#. For drawing a real number from a uniform distribution between 0 and 1, use ``Real Random::getUniformReal01()``.
#. For drawing a real number from a uniform distribution between -1 and 1, use ``Real Random::getUniformReal11()``.
#. For drawing a real number from a normal distribution centered at 0 and with a variance of 1, use ``Real Random::getNormal01()``.
#. For drawing an integer from a Poisson distribution with a specified mean, use ``T Random::getPoisson<T>(const Real a_mean)`` where ``T`` is an integer type.
#. For drawing a random direction in space, use ``RealVect Random::getDirection()``.
   The implementation uses the Marsaglia algorithm for drawing coordinates uniformly distributed over the unit sphere.


Setting the seed
----------------

By default, the random number generator is seeded with the MPI rank, which we do to avoid having the MPI ranks producing the same number sequences.
``Driver`` (see :ref:`Chap:Driver`) will seed the random number generator, and user can override the seed by setting ``Random.seed = <number>`` in the input script.
If the user sets ``<number> < 0`` then a random seed will be produced based on the elapsed CPU clock time.
If running with MPI, this seed is obtained by only one of the MPI ranks, and this seed is then broadcast to all the other ranks.
The other ranks will then increment the seed by their own MPI rank number so that each MPI rank gets a unique seed.

