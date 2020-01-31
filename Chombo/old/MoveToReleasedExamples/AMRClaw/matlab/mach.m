function mach = mach(data)
% MACH computes the mach number from conserved quantities.
%
%     MACH computes the mach number from Euler data in N dimension.
%     MACH assumes that the data  contains density in first column,
%     energy in last column, and components of momentum in columns
%     (2:DIM).
%
%      This function will be called from PLOTFRAME<N> if the user has set
%      the variable 'UserVariable' = 1 and 'UserVariableFile' = 'mach'.
%      These plotting parameters can be set in the SETPLOT<N> file.
%
%     gamma = 1.4 is hardwired here, but you can change this or modify
%     to read in the proper value from setprob.data, for example.
%
%  See also PRESSURE, XVELOCITY, YVELOCITY, ZVELOCITY, SETPLOT.

gamma = 1.4;
rho = data(:,1);
energy = data(:,end);
mom = data(:,2:end-1);
mom2 = mom .* mom;
kinetic = 0.5 * sum(mom2,2) ./ rho;
pressure = (gamma-1) * (energy - kinetic);
c2 = (gamma*pressure./rho);
speed2 = sum(mom2,2) ./ (rho.^2);
mach = sqrt(speed2 ./ c2);
