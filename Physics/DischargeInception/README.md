# Physics/DischargeInception
This physics module predicts discharge inception. 
The user must supply transport coefficients, geometry, and electron sources. 

## Setting up a new application
To set up a new problem, use the Python script. For example:

```shell
python3 setup.py -base_dir=MyApplications -app_name=MyDischargeInception -geometry=Cylinder
```

The application will then be installed to $DISCHARGE_HOME/MyApplications/MyDischargeInception.
