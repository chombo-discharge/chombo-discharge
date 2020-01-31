% Routines for plotting data output from Clawpack package
%
% Basic Clawpack plotting routines and options.
%
%      plotclaw1              - plot 1d results
%      plotclaw2              - plot 2d results
%      plotclaw3              - plot 3d results
%      setplot                - general help on setplot1, setplot2 and setplot3.
%
% Reading data, stepping through files, displaying and printing
%
%      readamrdata            - reads output files produced by Clawpack.
%      queryframe             - querys user for next frame information
%      printgif               - prints a gif file from current figure
%      printjpg               - prints a jpg file from current figure
%      setopengl              - Sets OpenGL renderer
%      setviews               - sets pre-defined viewing angles for 3d plots
%
% Defined functions for data analysis
%
%      pressure               - returns pressure given input data
%      xvelocity              - returns x velocity
%      yvelocity              - returns y velocity
%      zvelocity              - returns z velocity
%      mach                   - return mach number.
%
% General graph properties for 1d plots and line plots
%
%      setplotstyle           - sets symbols and colors for line/scatter plots
%      plotframe1ez           - functional form of plotframe1.
%      map1d                  - user defined function for creating 1d line plots
%      map1d_help             - help on creating 1d line plots from 2d/3d data
%      getlegendinfo          - returns legend information on line plots.
%
% Adding and modifying contour lines.
%
%      drawcontourlines       - draws contour lines on slices.
%      showcontourlines       - shows contour lines
%      hidecontourlines       - hides contour lines
%      setcontourlinecolor    - sets color of contour lines
%
% Modifying appearance of grid lines, meshes, patch borders (2d) and cubes (3d)
%
%      showgridlines          - shows computational grid
%      hidegridlines          - hides computational grid
%      setgridlinecolor       - sets color of grid lines
%
%      showpatchborders       - shows patch borders
%      hidepatchborders       - hides patch borders
%      setpatchbordercolor    - sets color of patch borders.
%
%      showmesh               - shows a coarsened mesh on specified levels
%      hidemesh               - hides a coarsened mesh on specified levels
%      setmeshcolor           - sets color of mesh.
%
%      showcubes              - shows 3d amr patch cube borders
%      hidecubes              - hides 3d amr patch cube borders
%      setcubecolor           - sets the color of 3d patch borders (cubes)
%
% Hiding/showing levels and slices for 2d/3d plots
%
%      showslices             - shows slices
%      hideslices             - hides slices
%      sliceloop              - loop over slices on 3d plots.
%
%      showlevels             - shows specified levels
%      hidelevels             - hides specified levels
%
% Modifying colors and alpha values of 2d/3d color slice plots
%
%      setcolors              - gives user control over how colors are set
%      setslicecolor          - sets color of slice/manifold
%      setslicealpha          - set transparency value of slice/manifold
%      yrbcolormap            - yellow/red/blue colormap
%      redwhite               - red/white colormap
%      rybcolormap            - red/yellow/blue colormap
%
% User defined functions for mapped grids
%
%      mappedgrids            - help on using mapped grids
%      mapc2p                 - User defined function for mapped grids
%      mapc2p_help            - help on using mapped grids
%      getblocknumber         - get block number for multi-block calculations
%
% User defined functions for 2d manifolds
%
%      manifold               - information on how to plot a manifold.
%      mapc2m                 - user-defined function for plotting manifolds
%      mapc2m_help            - help on plotting manifold data.
%      projectcontours        - projects contour lines to user-specified plane.
%
% Isosurface specific routines for 3d
%
%      showsurfs              - shows isosurfaces created with ISOSURFVALUES
%      hidesurfs              - hides isosurfaces created with ISOSURFVALUES
%      surfloop               - loop over isosurfaces
%
%      showsurflevels         - show isosurfaces at specified AMR levels
%      hidesurflevels         - hide isosurfaces at specified AMR levels
%      reducesurf             - reduces number of faces on isosurface.
%
%      showsurfmesh           - shows isosurface mesh
%      hidesurfmesh           - hides isosurface mesh
%      setsurfalpha           - sets isosurface transparency value
%      setsurfcolor           - sets isosurface color
%
% ChomboClaw related graphics
%
%      chomboclaw             - help for graphing using Chombo hdf5 output.
%      chombovis              - script for setting ChomboVis features.
%
% Type 'help' on any one of the individual topics above for more help.
%
