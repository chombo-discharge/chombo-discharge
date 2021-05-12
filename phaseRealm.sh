# # phase_realm begins here
# # Change inclusion guards. 
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
#     sed -i 's/\#include \"phase_realm.H\"/\#include \"CD_PhaseRealm.H\"/g' $i
#     sed -i 's/\#include <phase_realm.H>/\#include <CD_PhaseRealm.H>/g' $i
# done

# # Find and replace all instances of amr_mesh everywhere. This is now AmrMesh. 
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
#     sed -i 's/phase_realm/PhaseRealm/g' $i
# done

# # Move files
# mv src/AmrMesh/phase_realm.H    src/AmrMesh/CD_PhaseRealm.H
# mv src/AmrMesh/phase_realm.cpp  src/AmrMesh/CD_PhaseRealm.cpp


# # Update function signatures
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
#     sed -i 's/s_eb_NonConservativeDivergenceStencil/s_noncons_div/g' $i
#     sed -i 's/regrid_base/regridBase/g' $i
#     sed -i 's/query_operator/queryOperator/g' $i
# done

# # ----- realm begins here
# # Change inclusion guards. 
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
#     sed -i 's/\#include \"realm.H\"/\#include \"CD_Realm.H\"/g' $i
#     sed -i 's/\#include <realm.H>/\#include <CD_Realm.H>/g' $i
# done

# # Find and replace all instances of amr_mesh everywhere. This is now AmrMesh. 
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
#     sed -i 's/const static std::string primal/const static std::string Primal/g' $i
#     sed -i 's/realm/Realm/g' $i
#     sed -i 's/::primal/::Primal/g' $i
#     sed -i 's/\"primal\"/\"Primal\"/g' $i
# done

# # Move files
# mv src/AmrMesh/realm.H    src/AmrMesh/CD_Realm.H
# mv src/AmrMesh/realm.cpp  src/AmrMesh/CD_Realm.cpp

# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
#     sed -i 's/query_mask/queryMask/g' $i
# done


for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/get_realm/getRealm/g' $i
    sed -i 's/get_ebis/getEBIndexSpace/g' $i
    sed -i 's/get_gradsten/getGradientStencils/g' $i
done
