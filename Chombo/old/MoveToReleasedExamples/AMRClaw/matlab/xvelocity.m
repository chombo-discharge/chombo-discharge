function u = xvelocity(data)
%
% XVELOCITY computes the u-velocity from conserved quantities.
%
%      U = XVELOCITY(DATA) computes the u-velocity from Euler data
%      in N dimensions. XVELOCITY assumes data contains density in
%      first column, energy in last column, and components of momentum
%      in columns (1:DIM).
%
%      This function will be called from PLOTFRAME<N> if the user has set
%      the variable 'UserVariable' = 1 and 'UserVariableFile' = 'xvelocity'.
%      These plotting parameters can be set in the SETPLOT<N> file.
%
% See also PRESSURE, YVELOCITY, ZVELOCITY, MACH, SETPLOT.


u = data(:,2)./data(:,1);
