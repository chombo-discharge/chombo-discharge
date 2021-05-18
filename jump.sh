# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include <jump_bc.H>/\#include <CD_JumpBc.H>/g' $i
    sed -i 's/\#include \"jump_bc.H\"/\#include <CD_JumpBc.H>/g' $i
    sed -i 's/\#include <jump_bcI.H>/\#include <CD_JumpBcImplem.H>/g' $i
    sed -i 's/\#include \"jump_bcI.H\"/\#include <CD_JumpBcImplem.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/jump_bc/JumpBc/g' $i
    sed -i 's/match_bc/matchBc/g' $i
    sed -i 's/compute_dphidn/computeDnPhi/g' $i
    sed -i 's/get_stencils/getStencils/g' $i
    sed -i 's/get_avgStencils/getAvgStencils/g' $i
    sed -i 's/get_weights/getWeights/g' $i
    sed -i 's/get_avgWeights/getAvgWeights/g' $i
    sed -i 's/get_bco/getAvgBco/g' $i
    sed -i 's/get_inhomo/getInhomogeneous/g' $i
    sed -i 's/get_homog/getHomogeneous/g' $i
    sed -i 's/define_vofiter/defineVofIterator/g' $i
    sed -i 's/set_bco/setBcoefficient/g' $i
    sed -i 's/get_second_order_sten/getSecondOrderStencil/g' $i
    sed -i 's/get_first_order_sten/getFirstOrderStencil/g' $i
    sed -i 's/compute_avg_jump/computeAvgJump/g' $i
done

# Move files
mv Source/Elliptic/jump_bc.H     Source/Elliptic/CD_JumpBc.H
mv Source/Elliptic/jump_bcI.H    Source/Elliptic/CD_JumpBcImplem.H
mv Source/Elliptic/jump_bc.cpp   Source/Elliptic/CD_JumpBc.cpp

