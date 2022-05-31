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

## Building the Sphinx documentation locally
We use sphinx for building the user documentation. 
sphinx-autobuild is available on [PyPI](https://pypi.org/p/sphinx-autobuild/).
It can be installed using pip:

```
pip install sphinx-autobuild
```

### Manual build
To build the documentation locally, navigate to the Sphinx subfolder and run

```
make html
```

for the HTML documentation and

```
make latexpdf
```

for the LaTeX/PDF documentation.

Source files will end up in Sphinx/build/html and Sphinx/build/latex, respectively. 

### Adding Sphinx changes. 
To add changes to the Sphinx documentation, make the changes in the relevant Sphinx/source folder.
To build the documentation, either build it locally as above or put Sphinx in auto-build mode:

```
cd Sphinx
sphinx-autobuild source/ build/html
```

This will start a server at http://127.0.0.1:8000 and start watching for changes in the `source` directory.
When a change is detected, the documentation is rebuilt and any open browser windows are reloaded automatically. `KeyboardInterrupt` (<kbd>ctrl</kbd>+<kbd>c</kbd>) will stop the server.

