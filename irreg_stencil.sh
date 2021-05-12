# # Inclusion guards
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
#     sed -i 's/\#include \"irreg_stencil.H\"/\#include \"CD_IrregStencil.H\"/g' $i
#     sed -i 's/\#include <irreg_stencil.H>/\#include <CD_IrregStencil.H>/g' $i
# done

# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
#     sed -i 's/irreg_stencil/IrregStencil/g' $i
#     sed -i 's/stencil_type/IrregStencil::StencilType/g' $i
# done

# # Move files
# mv src/AmrMesh/irreg_stencil.H    src/AmrMesh/CD_IrregStencil.H
# mv src/AmrMesh/irreg_stencil.cpp  src/AmrMesh/CD_IrregStencil.cpp
    
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/build_stencil/buildStencil/g' $i
done
