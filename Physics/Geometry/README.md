# Physics/Geometry
This physics module only sets up a geometry -- it does include any solvers whatsoever and so all TimeStepper routines are empty.
This is typically used when developing/testing a new geometry. 

The source files consist of the following:

* **CD_GeometryStepper.H/cpp** Implementation of TimeStepper -- does not provide ANY solver functionality and can only instantiate a geometry. 

## Setting up a new application
To set up a new problem, use the Python script. For example:

```shell
python setup.py -base_dir=/home/foo/MyApplications -app_name=MyGeometry -geometry=Vessel
```

To install within chombo-discharge:

```shell
python setup.py -base_dir=$DISCHARGE_HOME/MyApplications -app_name=MyGeometry -geometry=Vessel
```

The application will be installed to $DISCHARGE_HOME/MyApplications/MyGeometry.

## Modifying the application
Users are free to modify this application.
