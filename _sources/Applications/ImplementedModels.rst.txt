.. _Chap:ImplementedModels:

Implemented models
==================

Various models that use the ``chombo-discharge`` solvers and functionality have been implemented and are shipped with ``chombo-discharge``.
These models implement ``TimeStepper`` for advancing various types of equations, and have been set up both for unit testing and as building blocks for more complex applications. 
The models reside in :file:`$DICHARGE_HOME/Physics` and are supplemented with Python tools for quickly setting up new applications that use the same code.
In general, users are encouraged to modify these models and tailor them for their own applications.
