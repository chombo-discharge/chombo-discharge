if (PlotType <= 3)
  yrbcolormap;
  grid on;
  axis([-3 3 -1 1 -1 1]);
  daspect([1 1 1]);
  caxis([1 1.25]);

  hideslices;
  showslices('x',2);
  showslices('y',2);
  showslices('z',2);
  showgridlines(1:2);

  showpatchborders;
  str = sprintf('Simple wave at t = %4.2f',t);
  title(str,'FontSize',16);
  xlabel('X','FontSize',16);
  ylabel('Y','FontSize',16);
  set(gca,'FontSize',14);
  hold off;
end;

if (exist('hcam'))
  if (ishandle(hcam))
    delete(hcam);
  end;
end;
hcam = camlight;

clear afterframe;
clear map1d;
