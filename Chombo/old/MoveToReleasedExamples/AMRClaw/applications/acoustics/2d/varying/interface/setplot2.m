
%  setplot2.m
%  called in plotclaw1 before plotting to set various parameters

OutputFlag = 'chombo';
plot_interval = 5;
plot_prefix = 'plotNEW';

PlotType = 1;                % type of plot to produce:
			     % 1 = pseudo-color (pcolor)
			     % 2 = contour
			     % 3 = Schlieren
			     % 4 = scatter plot of q vs. r

mq = 1;                      % which component of q to plot
UserVariable = 0;            % set to 1 to specify a user-defined variable
UserVariableFile = ' ';      % name of m-file mapping data to q
MappedGrid = 0;              % set to 1 if mapc2p.m exists for nonuniform grid
MaxFrames = 5000;            % max number of frames to loop over
MaxLevels = 6;
PlotData =  [1 1 1 1 1 1];   % Data on refinement level k is plotted only if
			     % k'th component is nonzero
PlotGrid =  [0 0 0 0 0 0];   % Plot grid lines on each level?

PlotGridEdges =  [1 1 1 1 1 1];  % Plot edges of patches of each grid at
                                 % this level?

%---------------------------------

% for contour plots (PlotType==2):
cv = linspace(-.1,0.1,15);    % contour levels
cv(cv == 0) = [];
ContourValues = cv;

%---------------------------------

LineStyle = setplotstyle('ro','gx');
UserMap1d = 1;
