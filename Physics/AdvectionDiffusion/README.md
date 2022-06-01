# Physics/AdvectionDiffusion
This physics module solves for an advection-diffusion process of a single scalar quantity. This module contains files for setting up the initial conditions
and advected species, basic integrators, and a cell tagger for refining grid cells. 

The source files consist of the following:

* **CD_AdvectionDiffusionSpecies.H** Implementation of CdrSpecies, for setting up initial conditions and turning on/off advection and diffusion.
* **CD_AdvectionDiffusionTagger.H/cpp**  Implementation of CellTagger, for flagging cells for refinement and coarsening.
* **CD_AdvectionDiffusionStepper.H/cpp** Implementation of TimeStepper, for advancing the advection-diffusion equation. 

# Setting up a new problem
To set up a new problem, use the Python script. For example:

```shell
./setup.py -base_dir=MyApplications -app_name=MyAdvectionDiffusion -geometry=CoaxialCable
```

The application will be installed to $DISCHARGE_HOME/MyApplications/MyAdvectionDiffusion.
The user will need to modify the geometry and set the initial conditions through the inputs file. 
