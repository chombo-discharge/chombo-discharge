global map_choice;
global skewness x_map;
global ubar vbar wbar;

fid = fopen('setprob.data','r');
gma = fscanf(fid,'%g %g %g',1);  s = fscanf(fid,'%s',1);
map_choice = fscanf(fid,'%d',1); s = fscanf(fid,'%s',1);
skewness = fscanf(fid,'%g',1); s = fscanf(fid,'%s',1);
x_map = fscanf(fid,'%g %g %g',3); s = fscanf(fid,'%s',1);
x_map = x_map/sqrt(dot(x_map,x_map));

disp([ubar vbar wbar]);
disp(map_choice);
disp(skewness);
disp(x_map);

fclose(fid);
