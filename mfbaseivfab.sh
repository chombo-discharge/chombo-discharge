#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"real_box.H\"/\#include <CD_RealBox.H>/g' $i
    sed -i 's/\#include <real_box.H>/\#include <CD_RealBox.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/real_box/RealBox/g' $i
    sed -i 's/get_lo/getLo/g' $i
    sed -i 's/get_hi/getHi/g' $i
    sed -i 's/get_corners/getCorners/g' $i
    sed -i 's/is_point_inside/isPointInside/g' $i
    sed -i 's/is_box_inside/isBoxInside/g' $i
done

# Move files
mv Source/Utilities/real_box.H      Source/Utilities/CD_RealBox.H
mv Source/Utilities/real_box.cpp    Source/Utilities/CD_RealBox.cpp
