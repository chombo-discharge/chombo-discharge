function pressure = pressure(data)
% PRESSURE computes the pressure from conserved quantities.
%
%      p = PRESSURE(data) computes the pressure from Euler data
%      in N dimensions.  PRESSURE assumes data contains density in
%      first column, energy in last column, and components of momentum
%      in columns (1:DIM).
%
%      This function will be called from PLOTFRAME<N> if the user has set
%      the variable 'UserVariable' = 1 and 'UserVariableFile' = 'pressure'.
%      These plotting parameters can be set in the SETPLOT<N> file.
%
%      gamma = 1.4 is hardwired here, but you can change this or modify
%      to read in the proper value from setprob.data, for example.
%
% See also XVELOCITY, YVELOCITY, ZVELOCITY, MACH, SETPLOT.
%

gamma = 1.4;
rho = data(:,1);
energy = data(:,end);
mom = data(:,2:end-1);
mom2 = mom .* mom;
kinetic = 0.5 * sum(mom2,2) ./ rho;
pressure = (gamma-1) * (energy - kinetic);
