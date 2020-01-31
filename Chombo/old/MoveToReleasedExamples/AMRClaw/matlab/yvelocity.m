function v = yvelocity(data)
%
% YVELOCITY computes the v-velocity from conserved quantities.
%
%      V = YVELOCITY(DATA) computes the v-velocity from Euler data
%      in N dimensions. YVELOCITY assumes data contains density in
%      first column, energy in last column, and components of momentum
%      in columns (1:DIM).
%
%      This function will be called from PLOTFRAME<N> if the user has set
%      the variable 'UserVariable' = 1 and 'UserVariableFile' = 'yvelocity'.
%      These plotting parameters can be set in the SETPLOT<N> file.
%
% See also PRESSURE, XVELOCITY, ZVELOCITY, MACH, SETPLOT.

v = data(:,3)./data(:,1);
