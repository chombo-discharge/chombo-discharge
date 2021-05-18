
for i in `find . -name "GNUmakefile" -type f`;
do
    sed -i 's/\/src/\/Source/g' $i
    sed -i 's/src\//Source\//g' $i
done

for i in `find . -name "*.py" -type f`;
do
    sed -i 's/\/src/\/Source/g' $i
    sed -i 's/src\//Source\//g' $i
done


mv src Source
