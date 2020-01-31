%     ChomboClaw is a set of C++ and Fortran routines freely available from
%     Lawrence Berkeley Laboratory for solving PDEs using an adaptive mesh
%     refinement framework.  The output produced from ChomboClaw
%     can now be plotted using the Clawpack Graphics routines.  Below, is a
%     description of graphics routines specific to Chombo.
%
%     Chombo and ChomboClaw produce output files using the HDF5 file format developed
%     the NCSA (see http://hdf.ncsa.uiuc.edu/HDF5).  The latest release of
%     Matlab, Version 6.5.1, release 13, has several routines for handling
%     HDF5 formatted files.  In the Clawpack graphics, we have added
%     routines which use these features to read the output from ChomboClaw
%     and produce data structures usuable by the Claw Graphics commands.
%
%     To plot output from a ChomboClaw run, there are only three parameters
%     that the user will need to set.  These paramters can be set in the
%     SETPLOT2 or SETPLOT3 files.  These parameters are :
%
%     OutputFlag        - set to 'chombo' for ChomboClaw output
%
%     plot_prefix       - set to prefix of file name of ChomboClaw
%                         output. For example, if a typical ChomboClaw file is
%                         called plotNEW0010.2d.hdf5, then the prefix is
%                         OutputPrefix = 'plotNEW'. By default, this prefix
%                         is 'pltstate'.   See also the input file variable
%                         'claw.plot_prefix' in the ChomboClaw User's Guide.
%
%     plot_interval     - set this to an integer specifying the plot
%                         interval.  For example, if output files are
%                         created every tenth time step, the plot_interval
%                         would be set to 10.  See also 'claw.plot_interval'
%                         in the ChomboClaw User's Guide.
%
%     chombovis         - this is a script which can be called from the top
%                         your setplot2 or setplot3 files to produce
%                         Matlab output that has the same "look and feel" as
%                         ChomboVis, the visualization package that comes
%                         with Chombo.
%
%     For more information on Chombo, see
%
%         http://seesar.lbl.gov/anag/chombo/index.html
%
%     For more information on ChomboClaw, contact Donna Calhoun at
%     calhoun@amath.washington.edu
%
%    See also SETPLOT, SETPLOT2, SETPLOT3, PLOTCLAW2, PLOTCLAW3, CHOMBOVIS.
%
