# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"irreg_amr_stencil.H\"/\#include \"CD_IrregAmrStencil.H\"/g' $i
    sed -i 's/\#include \"irreg_amr_stencilImplem.H\"/\#include \"CD_IrregAmrStencilImplem.H\"/g' $i
    
    sed -i 's/\#include <irreg_amr_stencil.H>/\#include <CD_IrregAmrStencil.H>/g' $i
    sed -i 's/\#include <irreg_amr_stencilImplem.H>/\#include <CD_IrregAmrStencilImplem.H>/g' $i
done

# Find and replace all instances of irreg_amr_stencil. This is now EbCentroidInterpolation
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/irreg_amr_stencil/IrregAmrStencil/g' $i
done

# Move files
mv src/AmrMesh/irreg_amr_stencil.H        src/AmrMesh/CD_IrregAmrStencil.H
mv src/AmrMesh/irreg_amr_stencilImplem.H  src/AmrMesh/CD_IrregAmrStencilImplem.H
    
