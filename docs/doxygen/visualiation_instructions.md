Visualization {#visualization-instructions}
=============

The code writes its output files to HDF5. These files are portable, binary data files which require specialized
visualization programs to plot. The recommended tool for this is to use VisIt visualization, which can be downloaded [here](https://wci.llnl.gov/simulation/computer-codes/visit/executables).

Three-dimensional data sets can become atrociously large. We are no strangers to data dumps in the range of tens of GB per time step; and local visualization of these are generally prohibitive.

VisIt can be run in server mode by allowing the local installation of VisIt to outsource the plot generation to a server. See
* [Cluster visualization](@ref visualization-cluster)