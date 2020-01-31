function [x,q] = map1d(xgrid,ygrid,qgrid)

[m,n] = size(qgrid);

x = linspace(0,5,100)';
y = x*0 + 2.5;
q = interp2(xgrid,ygrid,qgrid,x,y,'nearest*');
