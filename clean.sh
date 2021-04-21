# Insert driver.write_restart_file into all input scripts
for i in `find . -name "*.inputs" -type f`;
do
    sed -i 's/plasmac/chombo-discharge/g' $i
done

for i in `find . -name "*.ini" -type f`;
do
    sed -i 's/PlasmaC/chombo-discharge/g' $i
done

for i in `find . -name "*.py" -type f`;
do
    sed -i 's/PLASMAC_HOME/DISCHARGE_HOME/g' $i
done

for i in `find . -name "GNUmakefile" -type f`;
do
    sed -i 's/PLASMAC_HOME/DISCHARGE_HOME/g' $i
done

for i in `find . -name "*.py" -type f`;
do
    sed -i 's/plasmac/chombo-discharge/g' $i
done
