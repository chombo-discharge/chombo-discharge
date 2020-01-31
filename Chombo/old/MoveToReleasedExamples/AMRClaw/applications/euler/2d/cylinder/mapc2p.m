function [xp,yp] = mapc2p(xc,yc)
%
% specifies the mapping to curvilinear coordinates -- should be consistent
% with mapc2p.f
%
%

theta = 2*pi*yc;
r = 16*xc + 1;

xp = r .* cos(theta);
yp = r .* sin(theta);
