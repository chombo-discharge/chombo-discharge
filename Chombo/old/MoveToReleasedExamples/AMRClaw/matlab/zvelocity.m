function w = zvelocity(data)

%
% ZVELOCITY computes the w-velocity from conserved quantities.
%
%      W = ZVELOCITY(DATA) computes the w-velocity from Euler data
%      in N dimensions. ZVELOCITY assumes data contains density in
%      first column, energy in last column, and components of momentum
%      in columns (1:DIM).
%
%      This function will be called from PLOTFRAME<N> if the user has set
%      the variable 'UserVariable' = 1 and 'UserVariableFile' = 'zvelocity'.
%      These plotting parameters can be set in the SETPLOT<N> file.
%
% See also PRESSURE, XVELOCITY, YVELOCITY, MACH, SETPLOT.


w = data(:,4)./data(:,1);
