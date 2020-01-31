yrbcolormap;
axis([0 1 0 1 0 1]);
daspect([1 1 1]);
hideslices;
showslices('x',1);
showslices('y',1);
showslices('z',1);

showgridlines(1:2);

cv = linspace(-0.05,0.05,20);
cv(cv == 0) = [];
drawcontourlines(cv);

caxis([0 0.1])

chombovis;

clear afterframe;
