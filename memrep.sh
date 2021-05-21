#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"memrep.H\"/\#include <CD_MemoryReport.H>/g' $i
    sed -i 's/\#include <memrep.H>/\#include <CD_MemoryReport.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/memrep/MemoryReport/g' $i
    sed -i 's/get_memory/getMemoryUsage/g' $i
    sed -i 's/getMaxMin_memory/getMaxMinMemoryUsage/g' $i
done

# Move files
mv Source/Utilities/memrep.H    Source/Utilities/CD_MemoryReport.H
mv Source/Utilities/memrep.cpp  Source/Utilities/CD_MemoryReport.cpp
