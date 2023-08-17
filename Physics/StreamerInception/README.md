# Physics/StreamerInception
This physics module solves the streamer inception integral over the volume.
The user must supply the alpha and the geometry, and provide the voltages over which we evaluate the inception integral. 

## Setting up a new application
To set up a new problem, use the Python script. For example:

```shell
python3 setup.py -base_dir=MyApplications -app_name=MyStreamerInception -geometry=Cylinder
```

The application will then be installed to $DISCHARGE_HOME/MyApplications/MyStreamerInception.
