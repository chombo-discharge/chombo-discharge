function setmeshcolor(c)

%  SETMESHCOLOR sets the color of coarsened meshes.
%
%   SETMESHCOLOR(C) sets the color of coarsened meshes to C.
%   C must be a valid colorspec, either a string ('k','b','g', etc)
%   or an RGB triple.
%
%   See SHOWMESH, HIDEMESH.

sdirs = {'x', 'y', 'z'};
for idir = 1:3,
  slices = get_slices(sdirs{idir});
  for n = 1:length(slices),
    slice = slices{n};
    for level = 1:length(slice),
      pvec = slice{level};
      for k = 1:length(pvec),
	p = pvec(k);
	set_blocknumber(k);  % In case we are doing a multiblock plot
	udata = get(p,'UserData');
	if (isempty(udata.mesh.xlines) | isempty(udata.mesh.ylines))
	  continue;
	else
	  mesh = udata.mesh;
	end;
	xlines = mesh.xlines;
	ylines = mesh.ylines;
	lx = length(xlines);
	ly = length(ylines);
	set(xlines,'Color',c);
	set(ylines,'Color',c);
      end;
    end;
  end;
end;
