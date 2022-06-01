# chombo-discharge documentation quickstart

This is the chombo-discharge documentation.
chombo-discharge is hosted at https://github.com/chombo-discharge/chombo-discharge.

## Doxygen documentation
The doxygen documentation is available at https://chombo-discharge.github.io/doxygen/html/index.html
To build the doxygen documentation locally, navigate to the chombo-discharge root folder and type

```
doxygen Doxyfile
```

This will install the doxygen documentation in docs/doxygen. 

## Building the user documentation locally
We use sphinx for building the user documentation. 
sphinx-autobuild is available on [PyPI](https://pypi.org/p/sphinx-autobuild/).
It can be installed using pip:

```
pip install sphinx-autobuild
```

To build the documentation, run

```
sphinx-autobuild source/ ../build/html
```

This will start a server at http://127.0.0.1:8000 and start watching for changes in the `source` directory.
When a change is detected in `docs/`, the documentation is rebuilt and any open browser windows are reloaded automatically. `KeyboardInterrupt` (<kbd>ctrl</kbd>+<kbd>c</kbd>) will stop the server.

## Adding Sphinx changes. 
Use the github issue tracker and pull request system for contributions to the documentation.
To add to the Sphinx documentation, make your changes in the sphinx/source folder and run 'make clean github'.
This will compile the sphinx HTML documentation and append the result to the /docs folder (from which the documentation is built). 
