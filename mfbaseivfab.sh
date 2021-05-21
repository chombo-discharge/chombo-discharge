#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"MFQuadCFInterp.H\"/\#include <CD_MFQuadCFInterp.H>/g' $i
    sed -i 's/\#include <MFQuadCFInterp.H>/\#include <CD_MFQuadCFInterp.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/getNWOEBQuadCFInterp_ptr/getEBQuadCFInterpPointer/g' $i
    sed -i 's/getNWOEBQuadCFInterp/getEBQuadCFInterp/g' $i
done

# Move files
mv Source/Utilities/MFQuadCFInterp.H      Source/Multifluid/CD_MFQuadCFInterp.H
mv Source/Utilities/MFQuadCFInterp.cpp    Source/Multifluid/CD_MFQuadCFInterp.cpp
