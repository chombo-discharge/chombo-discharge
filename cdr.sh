# #Inclusion guards
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
#     sed -i 's/\#include \"mfis.H\"/\#include <CD_MultiFluidIndexSpace.H>/g' $i
#     sed -i 's/\#include <mfis.H>/\#include <CD_MultiFluidIndexSpace.H>/g' $i
# done

# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
#     sed -i 's/\_mfis/\_multiFluidIndexSpace/g' $i
#     sed -i 's/mfis/MultiFluidIndexSpace/g' $i
# done

# # Move files
# mv Source/Geometry/mfis.H   Source/Geometry/CD_MultiFluidIndexSpace.H
# mv Source/Geometry/mfis.cpp Source/Geometry/CD_MultiFluidIndexSpace.cpp

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/\Source\/CdrSolver/\/Source\/ConvectionDiffusionReaction/g' $i
    sed -i 's/\Source\/FieldSolver/\/Source\/Electrostatics/g' $i
done

mv Source/CdrSolver   Source/ConvectionDiffusionReaction
mv Source/FieldSolver Source/Electrostatics
