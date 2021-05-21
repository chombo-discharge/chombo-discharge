#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"poly.H\"/\#include <CD_PolyUtils.H>/g' $i
    sed -i 's/\#include <poly.H>/\#include <CD_PolyUtils.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/poly::/PolyUtils::/g' $i
    sed -i 's/brent_root_finder/brentRootFinder/g' $i
done

# Move files
mv Source/Utilities/poly.H      Source/Utilities/CD_PolyUtils.H
mv Source/Utilities/poly.cpp    Source/Utilities/CD_PolyUtils.cpp
