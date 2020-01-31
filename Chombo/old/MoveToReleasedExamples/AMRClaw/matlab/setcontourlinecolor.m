function setcontourlinecolor(c)

%  SETCONTOURLINECOLOR sets colors for contour lines.
%
%  SETCONTOURLINECOLOR(C) sets the color of contour lines to color C.
%  C should be a valid colorspec value, either a string indicating a
%  color value ('b','r','g', etc.) or an RGB triple.
%
%  See also SHOWCONTOURLINES, HIDECONTOURLINES, DRAWCONTOURLINES.
%

sdirs = {'x','y','z'};
for idir = 1:length(sdirs),
  slices = get_slices(sdirs{idir});
  for n = 1:length(slices),
    slice = slices{n};
    for level = 1:length(slice),
      pvec = slice{level};
      for k = 1:length(pvec),
	p = pvec(k);
	udata = get(p,'UserData');
	clines = udata.contourLines;
	set(clines,'Color',c);
      end;
    end;
  end;
end;
