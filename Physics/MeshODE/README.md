# Physics/MeshODE
This physics module solves for an ODE problem of a single scalar quantity. 

The source files consist of the following:

* **CD_MeshODEStepper.H** Implementation of TimeStepper, for advancing the equations of motion.

# Setting up a new problem
To set up a new problem, use the Python script. For example:

```shell
python setup.py -base_dir=/home/foo/MyApplications -app_name=MyProblem -geometry=Vessel
```

To install within chombo-discharge:

```shell
python setup.py -base_dir=$DISCHARGE_HOME/MyApplications -app_name=MyProblem -geometry=Vessel
```

The application will then be installed to $DISCHARGE_HOME/MyApplications/MyProblem.
The user will need to modify the geometry and set the initial conditions through the inputs file. 
