for i in `find . -name "GNUmakefile" -type f`;
do
    sed -i '/export/d' $i
done
