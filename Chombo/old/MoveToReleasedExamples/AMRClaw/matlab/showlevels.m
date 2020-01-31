function showlevels(level)

% SHOWLEVELS shows slice data
%
%       SHOWLEVELS(level) shows slice data at amr levels specified in vector
% 	  LEVEL.
%
%       SHOWLEVELS, by itself, shows slice data on all visible slices
%       at all amr levels.
%
%       NOTE : If the user has hidden all slices, it may not be possible
%       to show them again with SHOWLEVELS.  To recover visibility of
%       levels in this case, use SHOWSLICES.
%
%       See also HIDELEVELS, SHOWSLICES.


sdirs = {'x','y','z'};

off_on = {'off','on'};

for idir = 1:3,
  sdir = sdirs{idir};
  slices = get_slices(sdir);
  vs = get_visible_slices(sdir);
  vis_slices = ismember(1:length(slices),vs);
  for n = 1:length(slices),
    vstr = off_on{vis_slices(n)+1};
    slice = slices{n};
    if (nargin == 0)
      level = 1:length(slice);
    end;
    for l = 1:length(level),
      pvec = slice{level(l)};
      for k = 1:length(pvec),
	set(pvec(k),'Tag','on');

	% By setting visibility to vstr, rather than just 'on',
	% we don't show levels on slices that are not visible.
	set_patch_visibility(pvec(k),vstr);
      end;
    end;
    mask_patches_all(slice);
  end;
end;
