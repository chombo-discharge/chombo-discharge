function [rv,qv] = map1d(xcm,ycm,zcm,qcm)

global level;

[mx,my,mz] = size(xcm);

[xpm,ypm,zpm] = mapc2p(xcm,ycm,zcm);

if (level == 1)
  % qcm(xcm < 0) = nan;
end;

rv = reshape(xpm,mx*my*mz,1);
qv = reshape(qcm,mx*my*mz,1);
