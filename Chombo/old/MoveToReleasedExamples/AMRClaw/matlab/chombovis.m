% CHOMBOVIS sets up plots to look similar to ChomboVis plots produced by LBL.
%
%       CHOMBOVIS is a Matlab script that contains a series of graphics
%       commands that will set your axis colors, and background
%       colors to those of ChomboVis, the Chombo visualization package
%       developed at the Lawrence Berkeley Laboratory.
%
%       The script CHOMBOVIS sets the background color to black, turns off
%       the axis, and fixes the  viewing angle so that the aspect ratio
%       used for 3d plots doesn't change upon rotation.
%
%       Numeric axes values and labels can be turned on by using the 'axis on'
%       command.   By default, the axis is set to 'off';
%
%       For more information on Chombo and ChomboVis, see
%
%         http://seesar.lbl.gov/anag/chombo/index.html
%
%
%       See also CHOMBOCLAW.
%


% This fixes the aspect ratio under 3d rotations.
set(gca,'CameraViewAngleMode','manual');
set(gca,'CameraViewAngle',7);

% Set  background color to black
set(gcf,'Color','k');

% Set axis colors to gray (so it shows up on the black background).
clr_gray = [0.8 0.8 0.8];

h = get(gcf,'Children');
for i = 1:length(h),
  if (strcmp(get(h(i),'Type'),'axes') == 1)
    set(h(i),'Color','k');
    set(h(i),'XColor',clr_gray);
    set(h(i),'YColor',clr_gray);
    set(h(i),'ZColor',clr_gray);
    set(h(i),'Box','on')
  end;
end;

% set(gca,'Color','k');  % Set axis background itself to black
% set(gca,'XColor',clr_gray);  % Set x,y,z coordinate axes to gray.
% set(gca,'YColor',clr_gray);
% set(gca,'ZColor',clr_gray);
% set(gca,'Box','on');   % plot box.

% Default value.   Use 'axis on' to show axis.
axis off;
