if (PlotType==1)
  yrbcolormap
  caxis([0 0.2])
  colorbar
  hold on;
  plot([.5 .5],[0 1],'w')
  hold off
end
if (PlotType == 4)
  axis([0 1 0 1]);
  hold on
  fid = fopen('time.out','a');
  fprintf(fid,'%24.16e\n',t);
  fclose(fid);
  [qref_data,tref] = readamrdata(1,Frame/plot_interval,'./qref/');
  if (abs(tref - t) > 1e-6)
    error('Reference data time and current time are not compatible');
  end;

  hold on;
  [qref,xref,p] = plotframe1ez(qref_data,mq,'b-');
  hold off;
end;

axis square

clear afterframe;
