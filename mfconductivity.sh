# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include <CD_MfConductivityOp.H>/\#include <CD_MfHelmholtzOp.H>/g' $i
    sed -i 's/\#include <CD_MfConductivityOpFactory.H>/\#include <CD_MfHelmholtzOpFactory.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/MfConductivityOp/MfHelmholtzOp/g' $i
    sed -i 's/set_jump/setJump/g' $i
    sed -i 's/define_jump_stuff/defineJump/g' $i
    sed -i 's/define_multigrid_stuff/defineDeeperMultigrid/g' $i
    sed -i 's/averageDown_amr/averageDownAmr/g' $i
    sed -i 's/averageDown_mg/averageDownMG/g' $i
    sed -i 's/coarsen_coefficients/coarsenCoefficients/g' $i
    sed -i 's/set_bottom_drop/setBottomDrop/g' $i
    sed -i 's/set_relax_type/setRelaxType/g' $i
    sed -i 's/set_max_box_size/setMaxBoxSize/g' $i
    sed -i 's/set_ebbc_order/setEbBcOrder/g' $i
done

# Move files
mv Source/Elliptic/CD_MfConductivityOp.H           Source/Elliptic/CD_MfHelmholtzOp.H
mv Source/Elliptic/CD_MfConductivityOp.cpp         Source/Elliptic/CD_MfHelmholtzOp.cpp
mv Source/Elliptic/CD_MfConductivityOpFactory.H    Source/Elliptic/CD_MfHelmholtzOpFactory.H
mv Source/Elliptic/CD_MfConductivityOpFactory.cpp  Source/Elliptic/CD_MfHelmholtzOpFactory.cpp

