# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"lookup_table.H\"/\#include \"CD_LookupTable.H\"/g' $i
    sed -i 's/\#include <lookup_table.H>/\#include <CD_LookupTable.H>/g' $i
    sed -i 's/\#include \"lookup_tableI.H\"/\#include \"CD_LookupTableImplem.H\"/g' $i
    sed -i 's/\#include <lookup_tableI.H>/\#include <CD_LookupTableImplem.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/lookup_table/LookupTable/g' $i
    sed -i 's/make_uniform/makeUniform/g' $i
    sed -i 's/scale_x/scaleX/g' $i
    sed -i 's/scale_y/scaleY/g' $i
    sed -i 's/add_entry/addEntry/g' $i
    sed -i 's/add_table/addTable/g' $i
    sed -i 's/dump_table/dumpTable/g' $i
    sed -i 's/get_entry/getEntry/g' $i
    sed -i 's/direct_lookup/directLookup/g' $i
    sed -i 's/swap_xy/swapXY/g' $i
    sed -i 's/m_num_entries/m_numEntries/g' $i
done

# Move files
mv Source/Utilities/lookup_table.cpp  Source/Utilities/CD_LookupTable.cpp
mv Source/Utilities/lookup_table.H    Source/Utilities/CD_LookupTable.H
mv Source/Utilities/lookup_tableI.H   Source/Utilities/CD_LookupTableImplem.H

