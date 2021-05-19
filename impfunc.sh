# Insert footer in every .cpp and .H file
for i in `find . -name "*GNUmakefile" -type f`;
do
    var=`grep -n "(DISCHARGE_HOME)/Source/Geometry" $i | tail -1 | cut -f1 -d":"`
    echo $var
    sed -i $var'a\\t$(DISCHARGE_HOME)/Source/Multifluid \\' $i
    sed -i $var'a\\t$(DISCHARGE_HOME)/Source/ImplicitFunctions \\' $i
done
