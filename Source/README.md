# Source
This folder contains all chombo-discharge source code, excluding application codes (such as geometries and time integrators). 
The folder is organized as follows:

* AmrMesh Contains code for AMR and EB cores, including dual-mesh functionality and fundamental EB+AMR functionality.
* CellTagger Contains code for cell tagging.
* ConvectionDiffusionReaction Contains code for solving convection-diffusion-reaction problems with EB and AMR.
* Driver Contains the Driver class which run chombo-discharge simulations, and the TimeStepper class which is used in application codes.
* Electrostatics Contains source code for electrostatic field solves.
* Elliptic Contains source code for linear solvers (e.g., Helmholtz solvers).
* Geometry Contains source code for geometries. Includes declaration of an abstract geometry class (which the user can use in application codes), as well as support for polygon surfaces.
* ImplicitFunctions Contains source for for various useful implicit and signed distance functions.
* ItoDiffusion Contains source code for implementing drift-diffusion Brownian walkers (i.e., a particle code).
* Multifluid Contains source code for multifluid functionality.
* Particle Contains source code for using particles with EBs and AMR.
* RadiativeTransfer Contains source code for various radiative transfer solvers.
* SigmaSolver Contains source code for a surface charge solver.
* TracerParticles Contains source code for the tracer particle solver. 
* Utilities Contains various source code for useful utilities used in chombo-discharge source and application codes. 