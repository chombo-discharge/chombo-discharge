# Source
This folder contains all chombo-discharge source code, excluding application codes (such as geometries and time integrators). 
The folder is organized as follows:

* AmrMesh Contains code for AMR and EB cores, including dual-mesh functionality and fundamental EB+AMR functionality.
* CellTagger Contains code for cell tagging.
* ConvectionDiffusionReaction Contains code for solving convection-diffusion-reaction problems with EB and AMR.
* Driver Contains the Driver class which run chombo-discharge simulations, and the TimeStepper class which is used in application codes. 