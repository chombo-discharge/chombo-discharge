function mask_patches_all(slice)

% Internal Clawpack Graphics routine.

% This routine goes through and masks all patches,
% once the slice has been setup.

% First, lets figure out which levels are visible.

visible_levels = get_visible_levels(slice);

% Now, level visible_levels(l+1) should mask vis_level(l)

for l = 1:length(visible_levels),
  level = visible_levels(l);
  pvecl = slice{level};
  reset_level(pvecl);
  if (l < length(visible_levels))
    pvecu = slice{visible_levels(l+1)};
    if (~isempty(pvecu))
      for j = 1:length(pvecu),
	pj = pvecu(j);
	udata = get(pj,'UserData');
	xlow = udata.xmin;
	ylow = udata.ymin;
	zlow = udata.zmin;
	xhigh = udata.xmax;
	yhigh = udata.ymax;
	zhigh = udata.zmax;
	sdir = udata.sdir;

	for k = 1:length(pvecl)
	  pk = pvecl(k);
	  mask_patch(pk,sdir,xlow,xhigh,ylow,yhigh,zlow,zhigh);
	  mask_clines(pk,sdir,xlow,xhigh,ylow,yhigh,zlow,zhigh);
	  % mask_mesh(pk,sdir,xlow,xhigh,ylow,yhigh,zlow,zhigh);
	end;
      end;  % End loop on patches that need to be masked.
    end;
  end;
end;
