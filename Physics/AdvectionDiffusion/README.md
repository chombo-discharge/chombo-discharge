# Physics/AdvectionDiffusion
This physics module solves for an advection-diffusion process of a single scalar quantity. This module contains files for setting up the initial conditions
and advected species, basic integrators, and a cell tagger for refining grid cells. The source files consist of the following:

* **CD_AdvectionDiffusionSpecies.H/cpp** Implementation of CdrSpecies, for setting up initial conditions and turning on/off advection and diffusion.
* **CD_AdvectionDiffusionTagger.H/cpp**  Implementation of CellTagger, for flagging cells for refinement and coarsening.
* **CD_AdvectionDiffusionStepper.H/cpp** Implementation of TimeStepper, for advancing the advection-diffusion equation. 

## Setting up a new application
To set up a new problem, use

## Modifying the application
