# Physics/ItoKMC
This physics module solves the Ito-KMC formulation for electric discharges.

The folders contain:

* **PlasmaModels** Implementations of various plasma models. 
* **TimeSteppers** Implementations of the various time integration algorithms.
* **CellTaggers** Implementations of various cell taggers.

## Setting up a new application
To set up a new problem, use the Python script. For example:

```shell
python setup.py -base_dir=/home/foo/MyApplications -app_name=MyPlasmaModel -geometry=Vessel -physics=ItoKMCJSON
```

To install within chombo-discharge:

```shell
python setup.py -base_dir=$DISCHARGE_HOME/MyApplications -app_name=MyPlasmaModel -geometry=Vessel -physics=ItoKMCJSON
```

The application will then be installed to $DISCHARGE_HOME/MyApplications/MyPlasmaApplication.
The user will need to modify the geometry and set the initial conditions through the inputs file. 

## Modifying the application
Users are free to modify this application, e.g. adding new initial conditions and flow fields.
