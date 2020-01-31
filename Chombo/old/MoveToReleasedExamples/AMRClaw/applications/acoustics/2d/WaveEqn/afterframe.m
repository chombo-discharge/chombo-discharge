rybcolormap;

if (PlotType == 1)
  axis([0 5 0 5]);
  daspect([1 1 1]);
  caxis([-1 1]);
  % showgridlines(1:2);
  drawcontourlines(linspace(-1,1,20));
  colorbar;
elseif (PlotType == 4)
  Frame_1d = Frame/plot_interval;
  if (UserMap1d == 1)
    axis([0 5 -1 1]);
  else
    axis([0 3.5 -1 1]);
  end;
  [amr1d,t1d] =readamrdata(1,Frame_1d,'./1drad');
  if (abs(t - t1d) > 1e-12)
    fprintf('Times are incompatible; 1d time = %8.4f; 2d time = %8.4f\n',...
	t1,t);
    p = [];
  else
    [q1d,x1d,p] = plotframe1ez(amr1d,1,'k-');
    if (UserMap1d == 1)
      xd = get(p,'XData');
      xd = xd + 2.5;
      set(p,'XData',xd);
    end;
    set(p,'LineWidth',1);
  end;
  [h_amr,labels_amr] = getlegendinfo;
  if (UserMap1d == 1)
    pos = 3;
  else
    pos = 1;
  end;
  legend([h_amr,p],{labels_amr{:} sprintf('1d soln. (mx = %d)',...
      length(q1d))},pos);
  end


clear afterframe;
