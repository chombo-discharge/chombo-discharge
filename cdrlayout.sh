#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"cdr_layout.H\"/\#include <CD_CdrLayout.H>/g' $i
    sed -i 's/\#include <cdr_layout.H>/\#include <CD_CdrLayout.H>/g' $i
    sed -i 's/\#include \"cdr_layoutI.H\"/\#include <CD_CdrLayoutImplem.H>/g' $i
    sed -i 's/\#include <cdr_layoutI.H>/\#include <CD_CdrLayoutImplem.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/cdr\_layout/CdrLayout/g' $i
    sed -i 's/cdr\_iterator/CdrIterator/g' $i
    sed -i 's/cdr\_factory/CdrFactory/g' $i
    sed -i 's/get\_velocities/getVelocities/g' $i
done

# Remove files
rm Source/ConvectionDiffusionReaction/cdr_layout.H
rm Source/ConvectionDiffusionReaction/cdr_layoutI.H
