# Physics/Electrostatics
This physics module solves the Poisson equation. It does not feature mesh refinement, but this is possible to implement by adding a CellTagger to the module. 
See https://chombo-discharge.github.io/ElectrostaticsModel.html for implementation details.

The source files consist of the following:

* **CD_FieldStepper.H/cpp** Implementation of TimeStepper, for solving the Poisson equation. 

## Setting up a new application
To set up a new problem, use the Python script. For example:

```shell
./setup.py -base_dir=MyApplications -app_name=MyElectrostatics -geometry=CoaxialCable
```

The application will be installed to $DISCHARGE_HOME/MyApplications/MyElectrostatics.
The user will need to modify the geometry and set the initial conditions through the inputs file. 

## Modifying the application
Users are free to modify this application, e.g. adding support for mesh refinement or setting up more complex boundary conditions. 
