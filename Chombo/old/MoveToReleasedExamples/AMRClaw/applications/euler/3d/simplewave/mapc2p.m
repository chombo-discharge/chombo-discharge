function [xp,yp,zp] = mapc2p(xc,yc,zc)

global OutputFlag;

% xp = 2*xc - 3;
% yp = 3*yc - 1;
% zp = 5*zc - 1;


% # Rectangular grid with rotated inner sphere of radius xlen.

% xd = xd/sqrt(dot(xd,xd));

xd = [0.5 0.5 0.5];
x = xd(1);
y = xd(2);
z = xd(3);

skewness = 0.5;
radius = 0.8;
th = pi*skewness;

if (strcmp(OutputFlag,'chombo'))
  xc = xc - 3;
  yc = yc - 1;
  zc = zc - 1;
end

rcmat = sqrt(xc.^2 + yc.^2 + zc.^2);
amat = (1 + cos(pi*rcmat/radius))/2;
amat(rcmat > radius) = 0;
cmat = cos(amat*th);
smat = sin(amat*th);


[m,n,p] = size(xc);
for i = 1:m,
  for j = 1:n,
    for k = 1:p,
      rc = rcmat(i,j,k);
      a = amat(i,j,k);
      c = cmat(i,j,k);
      s = smat(i,j,k);

      Rot = [x^2*(1-c) + c    x*y*(1-c) - z*s x*z*(1-c) + y*s;
	x*y*(1-c) + z*s  y^2*(1-c) + c   y*z*(1-c) - x*s;
	x*z*(1-c) - y*s  y*z*(1-c) + x*s z^2*(1-c) + c  ];

      for ii = 1:3,
	xv(ii) = Rot(ii,1)*xc(i,j,k) + Rot(ii,2)*yc(i,j,k) + ...
	    Rot(ii,3)*zc(i,j,k);
      end;
      xp(i,j,k) = xv(1);
      yp(i,j,k) = xv(2);
      zp(i,j,k) = xv(3);
    end;
  end;
end;
