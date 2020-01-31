function hidelevels(level)

% HIDELEVELS hides slice data
%
%    HIDELEVELS, by itself,  hides slice data at all amr levels.
%
%    HIDELEVELS(level) hides slice data at amr levels specified in vector
% 	  LEVEL.
%
%    If you hide all levels, the slice is effectively not visible, and
%    you will not be able to show levels with SHOWLEVELS, as showlevels
%    only works on visible slices.  To see recover the slices, use
%    SHOWSLICES.
%
%    NOTE : If the user has hidden all levels using HIDELEVELS, it
%    may not be possible to show them again with SHOWLEVELS.
%    To recover visibility of levels in this case, use SHOWSLICES.
%
%    See also SHOWLEVELS, SHOWSLICES.
%

sdir = {'x','y','z'};

for idir = 1:3,
  slices = get_slices(sdir{idir});
  for n = 1:length(slices),
    slice = slices{n};
    if (nargin == 0)
      level = 1:length(slice);
    end;
    for l = 1:length(level),
      pvec = slice{level(l)};
      for k = 1:length(pvec),
	set(pvec(k),'Tag','off');
	set_patch_visibility(pvec(k),'off');
      end;
    end;
    mask_patches_all(slice);
  end;
end;
