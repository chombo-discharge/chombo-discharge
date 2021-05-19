#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"dcel_algorithms.H\"/\#include <CD_DcelAlgorithms.H>/g' $i
    sed -i 's/\#include <dcel_algorithms.H>/\#include <CD_DcelAlgorithms.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/dcel_algorithms/DcelAlgorithms/g' $i
done

# Move files

mv Source/Geometry/dcel_algorithms.H   Source/Geometry/CD_DcelAlgorithms.H
