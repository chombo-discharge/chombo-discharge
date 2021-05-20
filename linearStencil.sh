# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"LinearStencil.H\"/\#include \"CD_LinearStencil.H\"/g' $i
    sed -i 's/\#include <LinearStencil.H>/\#include <CD_LinearStencil.H>/g' $i
done

# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
#     sed -i 's/LinearStencil/LinearStencil/g' $i
# done

# Move files
mv Source/Utilities/LinearStencil.cpp Source/Utilities/CD_LinearStencil.cpp
mv Source/Utilities/LinearStencil.H   Source/Utilities/CD_LinearStencil.H

