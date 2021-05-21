#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"mfalias.H\"/\#include <CD_MultifluidAlias.H>/g' $i
    sed -i 's/\#include <mfalias.H>/\#include <CD_MultifluidAlias.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/mfalias/MultifluidAlias/g' $i
    sed -i 's/mf_cell_alias_factory/MfCellAliasFactory/g' $i
    sed -i 's/mf_flux_alias_factory/MfFluxAliasFactory/g' $i
    sed -i 's/mf_baseivfab_alias_factory/MfIVAliasFactory/g' $i
    sed -i 's/getMaxMin_memory/getMaxMinMemoryUsage/g' $i
done

# Move files
mv Source/Utilities/mfalias.H    Source/Utilities/CD_MultifluidAlias.H
