# Physics/CdrPlasma
This physics module solves the minimal plasma model

The folders contains

* **PlasmaModels** Implementations of various plasma models. 
* **TimeSteppers** Implementations of the various time integration algorithms.
* **CellTaggers** Implementations of various cell taggers.
* **Data** Various data (e.g. transport data). 


## Setting up a new application
To set up a new problem, use the Python script. For example:

```shell
./setup.py -base_dir=MyApplications -app_name=MyPlasmaModel -geometry=CoaxialCable
```

The application will be installed to $DISCHARGE_HOME/MyApplications/MyPlasmaApplication
The user will need to modify the geometry and set the initial conditions through the inputs file. 

## Modifying the application
Users are free to modify this application, e.g. adding new initial conditions and flow fields.
