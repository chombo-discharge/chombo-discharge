.. _Chap:Realm:

Realm
=====

The ``realm`` class is a class for centralizing EBAMR-related grids and operators for a specific AMR grid.
For example, a ``realm`` consists of a set of grids (i.e. a ``Vector<DisjointBoxLayout>``) as well as *operators*, e.g. functionality for filling ghost cells or averaging down a solution from a fine level to a coarse level. 

The terminology *dual grid* is used when more than one ``realm`` is used in a simulation.
With dual grid the user/developer has chosen to solve the equations of motion over a different set of ``DisjointBoxLayout`` on each level.
This approach is very useful when particle solvers are involved since users can quickly generate an Eulerian set of grids and a set of grids for the Lagrangian particles, and the grids can be load balanced separately.

In general, users will not interact with ``realm`` directly.
Every ``realm`` is owned by ``amr_mesh``, and the user will only interact with realms through the public ``amr_mesh`` interface.

Internally, an instantiation of ``realm`` contains the grids and the geometric information (e.g. EBISLayout), as well as any operators that the user has seen fit to *register*.
If a solver needs an operator for, say, ghost cell interpolation, the solver needs to *register* that operator through the ``amr_mesh`` public interface.
Operator registration is a run-time procedure.
Once an operator has been registered, ``realm`` will define those operators during e.g. regrids.
Run-time abortions with error messages are issued if an AMR operator is called for, but has not been registered. 
