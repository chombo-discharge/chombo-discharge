# Physics/RtPhysics
This physics module solves for a radiative transfer problem. This module contains files for setting up the initial conditions and integrators. It does not feature adaptive mesh refinement. 
See https://chombo-discharge.github.io/RtPhysicsModel.html for implementation details.

The source files consist of the following:

* **CD_RtPhysicsSpecies.H/cpp** Implementation of RtSpecies, for setting up a radiative transfer equation. 
* **CD_RtPhysicsStepper.H/cpp** Implementation of TimeStepper, for advancing the radiative transfer equation, either stationary or transient. 

## Setting up a new application
To set up a new problem, use the Python script. For example:

```shell
./setup.py -base_dir=myApplications -app_name=myRtPhysics -dim=2 -geometry=CoaxialCable -RtSolver=EddingtonSP1
```

Note that the user can choose between either using discrete or continuum models with the -RtSolver flag.
Acceptable options are EddingtonSP1 or McPhoto. 

The application will be installed to $DISCHARGE_HOME/myApplications/myRtPhysics.
The user will need to modify the geometry and set the initial conditions through the inputs file. 

## Modifying the application
By default, this application specifies a Gaussian source for the photons. 
Users are free to modify this application, e.g. adding other initial conditions.
