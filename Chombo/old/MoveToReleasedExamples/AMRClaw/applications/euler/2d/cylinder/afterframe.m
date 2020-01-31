% load default.cmap;
% colormap(default(:,1:3));
colormap('default');

axis([-5 5. -5 5])

caxis([0 23]);

daspect([1 1 1]);
cv = linspace(0,23,20);
drawcontourlines(cv);

% showgridlines(1);
% hidepatchborders;
% chombovis;

clear afterframe;
