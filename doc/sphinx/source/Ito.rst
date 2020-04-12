.. _Chap:ItoDiffusion:

Ito diffusion
=============

The Ito diffusion model advances computational particles as Brownian walkers with drift:

.. math::
   d\mathbf{X}_i = \mathbf{v}_idt + \sqrt{2D_i}\mathbf{W}_i dt,

where :math:`\mathbf{X}_i` is the spatial position of a particle :math:`i`, :math:`\mathbf{v}_i` is the drift coefficient and :math:`D_i` is the diffusion coefficient *in the continuum limit*.
That is, both :math:`\mathbf{v}_i` and :math:`D_i` are the quantities that appear in :ref:`Chap:CDR`.
The vector term :math:`\mathbf{W}_i` is a Gaussian random field with a mean value of 0 and standard deviation of 1.

The Îto particle
----------------

The Îto particle is a computational particle class in `PlasmaC` which can be used together with the particle tools in `Chombo`.
The following data fields are implemented in the particle:

.. code-block:: c++
   
   RealVect m_position;
   RealVect m_velocity;
   Real m_mass;
   Real m_diffusion;



Iterating over particles
------------------------

Limitations
-----------
