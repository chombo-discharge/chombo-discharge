# Physics/BrownianWalker
This physics module solves for advection-diffusion process of a single scalar quantity. This module contains files for setting up the initial conditions
and species, integrators, and a cell tagger for refining grid cells. See https://chombo-discharge.github.io/BrownianWalkerModel.html for implementation details.

The source files consist of the following:

* **CD_BrownianWalkerSpecies.H/cpp** Implementation of ItoSpecies, for setting up initial conditions and turning on/off advection and diffusion.
* **CD_BrownianWalkerTagger.H/cpp**  Implementation of CellTagger, for flagging cells for refinement and coarsening.
* **CD_BrownianWalkerStepper.H/cpp** Implementation of TimeStepper, for advancing the advection-diffusion equation. 

## Setting up a new application
To set up a new problem, use the Python script. For example:

```shell
./setup.py -base_dir=myApplications -app_name=myBrownianWalker -dim=2 -geometry=CoaxialCable
```

The application will be installed to $DISCHARGE_HOME/myApplications/myBrownianWalker.
The user will need to modify the geometry and set the initial conditions through the inputs file. 

## Modifying the application
The application is simply set up to advect and diffuse a scalar in a rotating flow.
Users are free to modify this application, e.g. adding new initial conditions and flow fields. 
