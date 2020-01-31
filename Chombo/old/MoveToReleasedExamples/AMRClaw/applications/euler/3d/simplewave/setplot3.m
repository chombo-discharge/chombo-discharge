%  setplot3.m
%  called in plotclaw3 before plotting to set various parameters

global OutputFlag;
OutputFlag = 'chombo';
plot_interval = 2;
plot_prefix = 'simple';
setopengl;

setviews  % set viewpoints so that view(xSlice), for example, can be used.

PlotType = 1;                % type of plot to produce:
			     % 1 = pcolor on slices (with optional contours)
			     % 2 = contour lines in 3d on transparent slices
			     % 3 = Schlieren plot on slices
			     % 4 = scatter plot of q vs. r
			     % 5 = isosurface plot (at given levels)

mq = 1;                      % which component of q to plot
UserVariable = 0;            % set to 1 to specify a user-defined variable
UserVariableFile = 'mach';      % name of m-file mapping data to q
MappedGrid = 1;              % set to 1 if mapc2p.m exists for nonuniform grid
MaxFrames = 1000;            % max number of frames to loop over
MaxLevels = 6;               % max number of AMR levels

PlotData =  [1 1 1 1 0 0];       % Data on refinement level k is plotted only
			         % if kth component is nonzero
PlotGrid =  [0 0 0 0 0 0];       % Plot grid lines on each level?
PlotGridEdges =  [0 0 0 0 0 0];  % Plot edges of patches of each grid at
                                 % this level on slices?
PlotCubeEdges = [0 0 0 0 0 0];   % Plot edges of cube of refinement patch at
                                 % this level?


ContourValues = linspace(0.95,1.3,20);

xSliceCoords = linspace(-3,3,3);
ySliceCoords = linspace(-1,1,3);
zSliceCoords = linspace(-1,1,3);

if (strcmp(OutputFlag,'chombo') == 1)
  xSliceCoords = xSliceCoords + 3;
  ySliceCoords = ySliceCoords + 1;
  zSliceCoords = zSliceCoords + 1;
end;

UserMap1d = 1;
ScatterStyle = setplotstyle('b+','g.','m.');
