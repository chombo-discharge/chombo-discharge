# Physics/TracerParticle
This physics module solves for an tracer particle solver.
This module contains files for setting up the initial conditions and advected species, basic integrators, and a cell tagger for refining grid cells.

The source files consist of the following:

* **CD_TracerParticleStepper.H**  Implementation of TimeStepper, for advancing a tracer particle solver. 

## Setting up a new application
To set up a new problem, use the Python script. For example:

```shell
python setup.py -base_dir=/home/foo/MyApplications -app_name=MyTracerParticle -geometry=Vessel
```

To install within chombo-discharge:

```shell
python setup.py -base_dir=$DISCHARGE_HOME/MyApplications -app_name=MyTracerParticle -geometry=Vessel
```

The application will then be installed to $DISCHARGE_HOME/MyApplications/TracerParticle.
The user will need to modify the geometry and set the initial conditions through the inputs file. 

## Modifying the application
The application is simply set up to track particles in various flow fields. 
Users are free to modify this application, e.g. adding new initial conditions and flow fields. 
