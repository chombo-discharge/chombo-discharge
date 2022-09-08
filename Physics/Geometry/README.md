# Physics/Geometry
This physics module only sets up a geometry -- it does include any solvers whatsoever and so all TimeStepper routines are empty.
This is typically used when developing/testing a new geometry. 
See https://chombo-discharge.github.io/Geometry.html for implementation details.

The source files consist of the following:

* **CD_GeometryStepper.H/cpp** Implementation of TimeStepper -- does not provide ANY solver functionality and can only instantiate a geometry. 

## Setting up a new application
To set up a new problem, use the Python script. For example:

```shell
./setup.py -base_dir=MyApplications -app_name=myGeometry -dim=2 -geometry=CoaxialCable
```

The application will be installed to $DISCHARGE_HOME/MyApplications/myGeometry.

## Modifying the application
Users are free to modify this application.
