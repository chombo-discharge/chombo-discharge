#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"data_ops.H\"/\#include <CD_DataOps.H>/g' $i
    sed -i 's/\#include <data_ops.H>/\#include <CD_DataOps.H>/g' $i
    sed -i 's/\#include \"data_opsI.H\"/\#include <CD_DataOpsImplem.H>/g' $i
    sed -i 's/\#include <data_opsI.H>/\#include <CD_DataOpsImplem.H>/g' $i
    sed -i 's/\#include \"data_opsF_F.H\"/\#include <CD_DataOpsF_F.H>/g' $i
    sed -i 's/\#include <data_opsF_F.H>/\#include <CD_DataOpsF_F.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/data_ops/DataOps/g' $i
    sed -i 's/set_value/setValue/g' $i
    sed -i 's/average_cell_to_face_allcomps/averageCellToFaceAllComps/g' $i
    sed -i 's/average_cell_to_face/averageCellToFace/g' $i
    sed -i 's/average_face_to_cell/averageFaceToCell/g' $i
    sed -i 's/dot_prod/dotProduct/g' $i
    sed -i 's/divide_scalar/divideByScalar/g' $i
    sed -i 's/get_max_min_norm/getMaxMinNorm/g' $i
    sed -i 's/get_max_min/getMaxMin/g' $i
    sed -i 's/kappa_sum/kappaSum/g' $i
    sed -i 's/kappa_scale/kappaScale/g' $i
    sed -i 's/gen_laplacian/genLaplacian/g' $i
    sed -i 's/flash_error/flashError/g' $i
    sed -i 's/multiply_scalar/multiplyScalar/g' $i
    sed -i 's/set_covered_valued/setCoveredValue/g' $i
    sed -i 's/square_root/squareRoot/g' $i
    sed -i 's/vector_length/vectorLength/g' $i
    sed -i 's/compute_min_valid_box/computeMinValidBox/g' $i
    sed -i 's/all_corners_inside_eb/allCornersInsideEb/g' $i
    sed -i 's/shift_corners/shiftCorners/g' $i
    sed -i 's/filter_smooth/filterSmooth/g' $i
    sed -i 's/compute_particle_weights/computeParticleWeights/g' $i
done

for i in `find . -type f \( -iname \*F_F.H  \)`; do
    sed -i 's/averageCellToFace/average_cell_to_face/g' $i
    sed -i 's/flashError/flash_error/g' $i
    sed -i 's/squareRoot/square_root/g' $i
    sed -i 's/vectorLength/vector_length/g' $i
    sed -i 's/genLaplacian/gen_laplacian/g' $i
    sed -i 's/filterSmooth/filter_smooth/g' $i
    sed -i 's/dotProduct/dot_product/g' $i
    sed -i 's/dot_productuct/dot_product/g' $i
    sed -i 's/averageFaceToCell/average_face_to_cell/g' $i
done

# Move files
mv Source/Utilities/data_ops.H    Source/Utilities/CD_DataOps.H
mv Source/Utilities/data_opsI.H   Source/Utilities/CD_DataOpsImplem.H
mv Source/Utilities/data_ops.cpp  Source/Utilities/CD_DataOps.cpp
mv Source/Utilities/data_opsF_F.H Source/Utilities/CD_DataOpsF_F.H
mv Source/Utilities/data_opsF.ChF Source/Utilities/CD_DataOpsF.ChF
